"""
Unit tests for representative genome mapping feature.

Tests the representative genome workflow: from model fields and metadata
handling through to classification and CLI commands. Representative mapping
allows ANI matrices to be built from a subset of genomes (one per species)
while alignment databases include all genomes, reducing ANI computation
from O(n^2) to O(k^2) where k is the number of species.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.mapping import app
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.models.genomes import AccessionList, GenomeAccession


# =============================================================================
# Helpers
# =============================================================================


def _make_taxonomy(genus: str, species: str, family: str = "Testaceae") -> str:
    """Build a minimal GTDB taxonomy string."""
    return f"d__Bacteria;f__{family};g__{genus};s__{species}"


# =============================================================================
# 1. GenomeAccession.is_representative field
# =============================================================================


class TestGenomeAccessionIsRepresentative:
    """Tests for the is_representative field on GenomeAccession."""

    def test_default_is_false(self) -> None:
        """is_representative should default to False."""
        acc = GenomeAccession(
            accession="GCF_001",
            gtdb_taxonomy="d__Bacteria;g__TestGenus;s__Test sp",
            species="Test sp",
        )
        assert acc.is_representative is False

    def test_explicit_true(self) -> None:
        """is_representative can be set to True."""
        acc = GenomeAccession(
            accession="GCF_001",
            gtdb_taxonomy="d__Bacteria;g__TestGenus;s__Test sp",
            species="Test sp",
            is_representative=True,
        )
        assert acc.is_representative is True

    def test_frozen_cannot_mutate(self) -> None:
        """is_representative should be immutable (frozen model)."""
        acc = GenomeAccession(
            accession="GCF_001",
            gtdb_taxonomy="d__Bacteria",
            species="Test sp",
            is_representative=False,
        )
        with pytest.raises(Exception):
            acc.is_representative = True


# =============================================================================
# 2. AccessionList.representative_map and to_metadata_tsv
# =============================================================================


class TestAccessionListRepresentativeMap:
    """Tests for AccessionList.representative_map and to_metadata_tsv output."""

    @pytest.fixture
    def accession_list_with_reps(self) -> AccessionList:
        """AccessionList with three species and representative assignments."""
        accessions = [
            GenomeAccession(
                accession="GCF_A1",
                gtdb_taxonomy=_make_taxonomy("GenusA", "Species alpha"),
                species="Species alpha",
                is_representative=True,
            ),
            GenomeAccession(
                accession="GCF_A2",
                gtdb_taxonomy=_make_taxonomy("GenusA", "Species alpha"),
                species="Species alpha",
                is_representative=False,
            ),
            GenomeAccession(
                accession="GCF_B1",
                gtdb_taxonomy=_make_taxonomy("GenusB", "Species beta"),
                species="Species beta",
                is_representative=True,
            ),
        ]
        return AccessionList(
            taxon="f__Testaceae",
            accessions=accessions,
            genus_counts={"GenusA": 2, "GenusB": 1},
            representative_map={
                "Species alpha": "GCF_A1",
                "Species beta": "GCF_B1",
            },
        )

    def test_representative_map_populated(self, accession_list_with_reps: AccessionList) -> None:
        """representative_map should contain species-to-accession entries."""
        rmap = accession_list_with_reps.representative_map
        assert rmap["Species alpha"] == "GCF_A1"
        assert rmap["Species beta"] == "GCF_B1"

    def test_representative_map_default_empty(self) -> None:
        """representative_map should default to an empty dict."""
        acc_list = AccessionList(taxon="test")
        assert acc_list.representative_map == {}

    def test_to_metadata_tsv_includes_representative_column(
        self,
        accession_list_with_reps: AccessionList,
        tmp_path: Path,
    ) -> None:
        """to_metadata_tsv should write a representative column."""
        out = tmp_path / "metadata.tsv"
        accession_list_with_reps.to_metadata_tsv(out)

        df = pl.read_csv(out, separator="\t")
        assert "representative" in df.columns

    def test_to_metadata_tsv_representative_values(
        self,
        accession_list_with_reps: AccessionList,
        tmp_path: Path,
    ) -> None:
        """Representative column should map each genome to its species rep."""
        out = tmp_path / "metadata.tsv"
        accession_list_with_reps.to_metadata_tsv(out)

        df = pl.read_csv(out, separator="\t")
        # GCF_A1 is representative of Species alpha - maps to itself
        row_a1 = df.filter(pl.col("accession") == "GCF_A1")
        assert row_a1["representative"][0] == "GCF_A1"

        # GCF_A2 is also Species alpha - maps to GCF_A1
        row_a2 = df.filter(pl.col("accession") == "GCF_A2")
        assert row_a2["representative"][0] == "GCF_A1"

        # GCF_B1 is representative of Species beta
        row_b1 = df.filter(pl.col("accession") == "GCF_B1")
        assert row_b1["representative"][0] == "GCF_B1"

    def test_to_metadata_tsv_no_representative_map(self, tmp_path: Path) -> None:
        """When representative_map is empty, each genome maps to itself."""
        accessions = [
            GenomeAccession(
                accession="GCF_X1",
                gtdb_taxonomy=_make_taxonomy("G", "S"),
                species="S",
            ),
        ]
        acc_list = AccessionList(
            taxon="test",
            accessions=accessions,
            representative_map={},
        )
        out = tmp_path / "metadata.tsv"
        acc_list.to_metadata_tsv(out)

        df = pl.read_csv(out, separator="\t")
        assert df["representative"][0] == "GCF_X1"


# =============================================================================
# 3. AccessionList.from_tsv with representative column
# =============================================================================


class TestAccessionListFromTsvRepresentative:
    """Tests for reading representative information from TSV files."""

    def test_from_tsv_reads_representative_column(self, tmp_path: Path) -> None:
        """from_tsv should populate representative_map from file."""
        tsv_content = (
            "accession\tgtdb_taxonomy\tspecies\trepresentative\n"
            "GCF_001\td__Bacteria;g__G;s__S1\tS1\tGCF_001\n"
            "GCF_002\td__Bacteria;g__G;s__S1\tS1\tGCF_001\n"
            "GCF_003\td__Bacteria;g__G;s__S2\tS2\tGCF_003\n"
        )
        tsv_path = tmp_path / "accessions.tsv"
        tsv_path.write_text(tsv_content)

        acc_list = AccessionList.from_tsv(tsv_path)

        assert acc_list.representative_map["S1"] == "GCF_001"
        assert acc_list.representative_map["S2"] == "GCF_003"

    def test_from_tsv_sets_is_representative(self, tmp_path: Path) -> None:
        """from_tsv should set is_representative=True for self-referencing rows."""
        tsv_content = (
            "accession\tgtdb_taxonomy\tspecies\trepresentative\n"
            "GCF_001\td__Bacteria\tS1\tGCF_001\n"
            "GCF_002\td__Bacteria\tS1\tGCF_001\n"
        )
        tsv_path = tmp_path / "accessions.tsv"
        tsv_path.write_text(tsv_content)

        acc_list = AccessionList.from_tsv(tsv_path)

        rep_flags = {a.accession: a.is_representative for a in acc_list.accessions}
        assert rep_flags["GCF_001"] is True
        assert rep_flags["GCF_002"] is False

    def test_from_tsv_no_representative_column(self, tmp_path: Path) -> None:
        """from_tsv should handle files without representative column."""
        tsv_content = (
            "accession\tgtdb_taxonomy\tspecies\n"
            "GCF_001\td__Bacteria\tS1\n"
        )
        tsv_path = tmp_path / "accessions.tsv"
        tsv_path.write_text(tsv_content)

        acc_list = AccessionList.from_tsv(tsv_path)

        assert acc_list.representative_map == {}
        assert acc_list.accessions[0].is_representative is False

    def test_roundtrip_metadata_preserves_representatives(self, tmp_path: Path) -> None:
        """Writing and reading metadata should preserve representative mapping."""
        accessions = [
            GenomeAccession(
                accession="GCF_R1",
                gtdb_taxonomy=_make_taxonomy("G1", "Sp1"),
                species="Sp1",
                is_representative=True,
            ),
            GenomeAccession(
                accession="GCF_R2",
                gtdb_taxonomy=_make_taxonomy("G1", "Sp1"),
                species="Sp1",
            ),
        ]
        original = AccessionList(
            taxon="test",
            accessions=accessions,
            representative_map={"Sp1": "GCF_R1"},
        )

        metadata_path = tmp_path / "metadata.tsv"
        original.to_metadata_tsv(metadata_path)

        loaded = AccessionList.from_tsv(metadata_path)
        assert loaded.representative_map["Sp1"] == "GCF_R1"


# =============================================================================
# 4. GenomeMetadata representative methods
# =============================================================================


class TestGenomeMetadataRepresentativeMethods:
    """Tests for GenomeMetadata representative-related methods."""

    @pytest.fixture
    def metadata_with_reps(self) -> GenomeMetadata:
        """GenomeMetadata with representative column."""
        df = pl.DataFrame({
            "accession": ["GCF_001", "GCF_002", "GCF_003", "GCF_004"],
            "species": ["SpA", "SpA", "SpB", "SpB"],
            "genus": ["GA", "GA", "GB", "GB"],
            "family": ["FA", "FA", "FA", "FA"],
            "representative": ["GCF_001", "GCF_001", "GCF_003", "GCF_003"],
        })
        return GenomeMetadata(df)

    @pytest.fixture
    def metadata_without_reps(self) -> GenomeMetadata:
        """GenomeMetadata without representative column."""
        df = pl.DataFrame({
            "accession": ["GCF_001", "GCF_002"],
            "species": ["SpA", "SpB"],
            "genus": ["GA", "GB"],
        })
        return GenomeMetadata(df)

    # -- get_representative --

    def test_get_representative_returns_rep(self, metadata_with_reps: GenomeMetadata) -> None:
        """get_representative should return the species representative."""
        assert metadata_with_reps.get_representative("GCF_002") == "GCF_001"

    def test_get_representative_self_for_rep(self, metadata_with_reps: GenomeMetadata) -> None:
        """get_representative should return self for a representative genome."""
        assert metadata_with_reps.get_representative("GCF_001") == "GCF_001"

    def test_get_representative_unknown_accession(
        self, metadata_with_reps: GenomeMetadata
    ) -> None:
        """get_representative should return accession itself for unknown genomes."""
        assert metadata_with_reps.get_representative("GCF_UNKNOWN") == "GCF_UNKNOWN"

    def test_get_representative_no_column(
        self, metadata_without_reps: GenomeMetadata
    ) -> None:
        """get_representative should return identity when no representative column."""
        assert metadata_without_reps.get_representative("GCF_001") == "GCF_001"

    # -- build_representative_mapping --

    def test_build_representative_mapping(self, metadata_with_reps: GenomeMetadata) -> None:
        """build_representative_mapping should return complete mapping dict."""
        mapping = metadata_with_reps.build_representative_mapping()
        assert mapping["GCF_001"] == "GCF_001"
        assert mapping["GCF_002"] == "GCF_001"
        assert mapping["GCF_003"] == "GCF_003"
        assert mapping["GCF_004"] == "GCF_003"

    def test_build_representative_mapping_identity(
        self, metadata_without_reps: GenomeMetadata
    ) -> None:
        """Without representative column, mapping should be identity."""
        mapping = metadata_without_reps.build_representative_mapping()
        assert mapping["GCF_001"] == "GCF_001"
        assert mapping["GCF_002"] == "GCF_002"

    # -- has_representatives --

    def test_has_representatives_true(self, metadata_with_reps: GenomeMetadata) -> None:
        """has_representatives should be True when non-self mappings exist."""
        assert metadata_with_reps.has_representatives is True

    def test_has_representatives_false_no_column(
        self, metadata_without_reps: GenomeMetadata
    ) -> None:
        """has_representatives should be False without representative column."""
        assert metadata_without_reps.has_representatives is False

    def test_has_representatives_false_all_self(self) -> None:
        """has_representatives should be False if all genomes map to themselves."""
        df = pl.DataFrame({
            "accession": ["GCF_001", "GCF_002"],
            "species": ["SpA", "SpB"],
            "genus": ["GA", "GB"],
            "representative": ["GCF_001", "GCF_002"],
        })
        meta = GenomeMetadata(df)
        assert meta.has_representatives is False

    # -- representative_count --

    def test_representative_count_with_reps(
        self, metadata_with_reps: GenomeMetadata
    ) -> None:
        """representative_count should return the number of unique representatives."""
        # 4 genomes, 2 representatives (GCF_001, GCF_003)
        assert metadata_with_reps.representative_count == 2

    def test_representative_count_no_column(
        self, metadata_without_reps: GenomeMetadata
    ) -> None:
        """representative_count should equal genome_count without representative column."""
        assert metadata_without_reps.representative_count == 2

    # -- Backwards compatibility --

    def test_backwards_compatibility_lookups(
        self, metadata_without_reps: GenomeMetadata
    ) -> None:
        """All representative methods should work when representative column is absent."""
        meta = metadata_without_reps
        # get_representative returns identity
        assert meta.get_representative("GCF_001") == "GCF_001"
        # build_representative_mapping returns identity
        mapping = meta.build_representative_mapping()
        assert all(k == v for k, v in mapping.items())
        # has_representatives is False
        assert meta.has_representatives is False
        # representative_count equals genome_count
        assert meta.representative_count == meta.genome_count


# =============================================================================
# 5. VectorizedClassifier with representative_mapping
# =============================================================================


class TestVectorizedClassifierWithRepresentativeMapping:
    """Tests for VectorizedClassifier using representative_mapping.

    Scenario: 10 genomes belonging to 3 species, with 3 representatives.
    The ANI matrix contains only the 3 representatives. BLAST output
    hits all 10 genomes. Classification should use representative ANI
    values while preserving actual hit genomes in best_match_genome.
    """

    @pytest.fixture
    def species_setup(self) -> dict:
        """Define species, genomes, and representative assignments.

        Returns a dictionary with species layout, representative mapping,
        ANI matrix, and BLAST DataFrame.
        """
        # Species layout: 3 species with varying numbers of genomes
        species_genomes = {
            "Species_alpha": ["GCF_A01", "GCF_A02", "GCF_A03", "GCF_A04"],
            "Species_beta": ["GCF_B01", "GCF_B02", "GCF_B03"],
            "Species_gamma": ["GCF_G01", "GCF_G02", "GCF_G03"],
        }

        # Representatives: one per species
        representatives = {
            "Species_alpha": "GCF_A01",
            "Species_beta": "GCF_B01",
            "Species_gamma": "GCF_G01",
        }

        # Build genome -> representative mapping (for all 10 genomes)
        representative_mapping: dict[str, str] = {}
        for species, genomes in species_genomes.items():
            rep = representatives[species]
            for g in genomes:
                representative_mapping[g] = rep

        # ANI matrix with only the 3 representatives
        ani_dict = {
            "GCF_A01": {"GCF_B01": 85.0, "GCF_G01": 75.0},
            "GCF_B01": {"GCF_A01": 85.0, "GCF_G01": 78.0},
            "GCF_G01": {"GCF_A01": 75.0, "GCF_B01": 78.0},
        }

        # BLAST DataFrame: reads hitting various genomes across all species.
        # Use pipe-separated sseqid format: "{accession}|{contig_id}"
        # Bitscores are kept close (within 5% of max) so that hits from
        # different species are above the default 95% bitscore threshold
        # and thus included in the ambiguous hit set.
        blast_rows = []
        all_genomes = []
        for genomes in species_genomes.values():
            all_genomes.extend(genomes)

        # Create 5 reads, each hitting all 10 genomes with similar bitscores
        for read_idx in range(5):
            for g_idx, genome in enumerate(all_genomes):
                # All species get similar bitscores (within 5% of 300)
                # to ensure cross-species ambiguity
                if genome.startswith("GCF_A"):
                    pident = 99.0 - g_idx * 0.1
                    bitscore = 300.0 - g_idx * 1.0
                elif genome.startswith("GCF_B"):
                    pident = 95.0 - (g_idx - 4) * 0.1
                    bitscore = 296.0 - (g_idx - 4) * 1.0
                else:
                    pident = 90.0 - (g_idx - 7) * 0.1
                    bitscore = 292.0 - (g_idx - 7) * 1.0

                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}|contig_1",
                    "pident": pident,
                    "length": 250,
                    "mismatch": int(250 * (1 - pident / 100)),
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 250,
                    "sstart": 1000,
                    "send": 1250,
                    "evalue": 1e-50,
                    "bitscore": bitscore,
                    "genome_name": genome,
                })

        blast_df = pl.DataFrame(blast_rows)

        return {
            "species_genomes": species_genomes,
            "representatives": representatives,
            "representative_mapping": representative_mapping,
            "ani_dict": ani_dict,
            "blast_df": blast_df,
        }

    def test_classification_with_representative_mapping(self, species_setup: dict) -> None:
        """Classification should succeed with representative mapping."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(species_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=species_setup["representative_mapping"],
        )

        result = classifier.classify_file(species_setup["blast_df"])

        assert isinstance(result, pl.DataFrame)
        assert len(result) == 5  # 5 reads

    def test_best_match_genome_is_actual_genome(self, species_setup: dict) -> None:
        """best_match_genome should be the actual hit genome, not the representative."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(species_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=species_setup["representative_mapping"],
        )

        result = classifier.classify_file(species_setup["blast_df"])

        best_genomes = set(result["best_match_genome"].to_list())
        rep_genomes = {"GCF_A01", "GCF_B01", "GCF_G01"}

        # The best genome should come from one of the actual genomes in the
        # BLAST output. In this scenario, the highest scoring genome for each
        # read is from Species alpha (GCF_A01 has the highest bitscore).
        # It may be GCF_A01 (which is also a representative) but should not
        # be forced to a representative.
        for genome in best_genomes:
            # Must be one of the 10 actual genomes, not some unknown
            all_genomes = []
            for genomes_list in species_setup["species_genomes"].values():
                all_genomes.extend(genomes_list)
            assert genome in all_genomes

    def test_classification_uses_representative_ani(self, species_setup: dict) -> None:
        """ANI lookups should use representative genomes for uncertainty calculation.

        When a read hits genomes from multiple species, the placement uncertainty
        should be derived from ANI values between the representatives in the
        matrix, not from the actual genomes (which are absent from the matrix).
        """
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(species_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=species_setup["representative_mapping"],
        )

        result = classifier.classify_file(species_setup["blast_df"])

        # Reads hit genomes from all 3 species (all within bitscore threshold).
        # The classifier should use ANI values from the representative matrix
        # to compute placement uncertainty. Without the mapping, ANI lookups
        # for non-representative genomes would return 0 (missing), giving
        # max uncertainty of 100%. With the mapping, ANI values of 85, 75, 78
        # should be used, yielding much lower uncertainty.
        uncertainties = result["placement_uncertainty"].to_list()
        for u in uncertainties:
            # ANI between reps is 75-85%, so uncertainty = 100 - ANI = 15-25%.
            # Without mapping it would be 100% (ANI = 0 for missing genomes).
            assert u < 30.0, (
                f"Uncertainty {u} suggests representative ANI was not used"
            )

    def test_classification_without_mapping_has_higher_uncertainty(
        self, species_setup: dict
    ) -> None:
        """Without mapping, many genomes miss ANI lookups, increasing uncertainty."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(species_setup["ani_dict"])

        # With mapping
        with_mapping = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=species_setup["representative_mapping"],
        )
        result_with = with_mapping.classify_file(species_setup["blast_df"])

        # Without mapping
        without_mapping = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=None,
        )
        result_without = without_mapping.classify_file(species_setup["blast_df"])

        mean_u_with = result_with["placement_uncertainty"].mean()
        mean_u_without = result_without["placement_uncertainty"].mean()

        # Without mapping, most genomes are absent from the ANI matrix, so
        # the classifier cannot look up ANI values and defaults to high
        # uncertainty. With mapping, ANI is resolved through representatives.
        assert mean_u_with < mean_u_without, (
            f"Mapping should reduce uncertainty: {mean_u_with} vs {mean_u_without}"
        )

    def test_all_expected_columns_present(self, species_setup: dict) -> None:
        """Output should contain all standard classification columns."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(species_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=species_setup["representative_mapping"],
        )

        result = classifier.classify_file(species_setup["blast_df"])

        expected_cols = {
            "read_id",
            "best_match_genome",
            "top_hit_identity",
            "novelty_index",
            "placement_uncertainty",
            "genus_uncertainty",
            "ambiguity_scope",
            "num_ambiguous_hits",
            "taxonomic_call",
            "diversity_status",
            "is_novel",
            "confidence_score",
        }
        assert expected_cols.issubset(set(result.columns))


# =============================================================================
# 6. validate_ani_genome_coverage with representative_mapping
# =============================================================================


class TestValidateAniGenomeCoverageWithRepresentativeMapping:
    """Tests for validate_ani_genome_coverage with representative_mapping."""

    def _write_blast_file(
        self, path: Path, genomes: list[str]
    ) -> None:
        """Write a minimal BLAST file with hits to the given genomes."""
        rows = []
        for i, genome in enumerate(genomes):
            rows.append({
                "qseqid": f"read_{i:03d}",
                "sseqid": f"{genome}|contig_1",
                "pident": 98.0,
                "length": 150,
                "mismatch": 3,
                "gapopen": 0,
                "qstart": 1,
                "qend": 150,
                "sstart": 1000,
                "send": 1150,
                "evalue": 1e-50,
                "bitscore": 250.0,
            })
        df = pl.DataFrame(rows)
        df.write_csv(path, separator="\t", include_header=False)

    def test_coverage_with_representative_mapping(self, tmp_path: Path) -> None:
        """With mapping, coverage should reflect representative presence in ANI."""
        from metadarkmatter.cli.score import validate_ani_genome_coverage

        # ANI matrix has only representatives
        ani_dict = {
            "GCF_REP1": {"GCF_REP2": 90.0},
            "GCF_REP2": {"GCF_REP1": 90.0},
        }
        ani_matrix = ANIMatrix(ani_dict)

        # BLAST file has non-representative genomes
        blast_path = tmp_path / "blast.tsv"
        self._write_blast_file(blast_path, ["GCF_001", "GCF_002", "GCF_003"])

        # Mapping: all 3 genomes map to 2 representatives
        rep_mapping = {
            "GCF_001": "GCF_REP1",
            "GCF_002": "GCF_REP1",
            "GCF_003": "GCF_REP2",
        }

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix, representative_mapping=rep_mapping,
        )

        # After mapping, check_genomes = {GCF_REP1, GCF_REP2}, both in ANI
        assert matched == 2
        assert total == 2
        assert pct == 100.0
        assert len(missing) == 0

    def test_coverage_without_mapping_shows_missing(self, tmp_path: Path) -> None:
        """Without mapping, raw genome names are checked against ANI."""
        from metadarkmatter.cli.score import validate_ani_genome_coverage

        ani_dict = {
            "GCF_REP1": {"GCF_REP2": 90.0},
            "GCF_REP2": {"GCF_REP1": 90.0},
        }
        ani_matrix = ANIMatrix(ani_dict)

        blast_path = tmp_path / "blast.tsv"
        self._write_blast_file(blast_path, ["GCF_001", "GCF_002"])

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix, representative_mapping=None,
        )

        # Without mapping, GCF_001 and GCF_002 are not in ANI
        assert matched == 0
        assert total == 2
        assert pct == 0.0
        assert "GCF_001" in missing
        assert "GCF_002" in missing

    def test_partial_representative_coverage(self, tmp_path: Path) -> None:
        """Some representatives missing from ANI should be reported."""
        from metadarkmatter.cli.score import validate_ani_genome_coverage

        ani_dict = {
            "GCF_REP1": {"GCF_REP1": 100.0},
        }
        ani_matrix = ANIMatrix(ani_dict)

        blast_path = tmp_path / "blast.tsv"
        self._write_blast_file(blast_path, ["GCF_001", "GCF_002"])

        rep_mapping = {
            "GCF_001": "GCF_REP1",
            "GCF_002": "GCF_REP_MISSING",
        }

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix, representative_mapping=rep_mapping,
        )

        assert matched == 1
        assert total == 2
        assert pct == 50.0
        assert "GCF_REP_MISSING" in missing


# =============================================================================
# 7. select-representatives CLI command
# =============================================================================


class TestSelectRepresentativesCLI:
    """Tests for the select-representatives CLI command."""

    @pytest.fixture
    def runner(self) -> CliRunner:
        """Create a Typer CLI runner."""
        return CliRunner()

    def _write_metadata(
        self,
        path: Path,
        *,
        include_representative: bool = False,
    ) -> None:
        """Write a test metadata TSV file."""
        data: dict[str, list[str]] = {
            "accession": ["GCF_B01", "GCF_A02", "GCF_A01", "GCF_B02"],
            "species": ["SpB", "SpA", "SpA", "SpB"],
            "genus": ["GB", "GA", "GA", "GB"],
            "family": ["FA", "FA", "FA", "FA"],
        }
        if include_representative:
            data["representative"] = ["GCF_B01", "GCF_A02", "GCF_A02", "GCF_B01"]

        df = pl.DataFrame(data)
        df.write_csv(path, separator="\t")

    def test_basic_representative_assignment(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should assign first alphabetical accession as representative."""
        metadata_path = tmp_path / "metadata.tsv"
        self._write_metadata(metadata_path)
        output_path = tmp_path / "output.tsv"

        result = runner.invoke(
            app,
            [
                "select-representatives",
                "--metadata", str(metadata_path),
                "--output", str(output_path),
            ],
        )

        assert result.exit_code == 0
        assert output_path.exists()

        df = pl.read_csv(output_path, separator="\t")
        assert "representative" in df.columns

        # SpA: first alphabetically = GCF_A01
        spa_reps = df.filter(pl.col("species") == "SpA")["representative"].to_list()
        assert all(r == "GCF_A01" for r in spa_reps)

        # SpB: first alphabetically = GCF_B01
        spb_reps = df.filter(pl.col("species") == "SpB")["representative"].to_list()
        assert all(r == "GCF_B01" for r in spb_reps)

    def test_preserves_existing_representatives(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should preserve self-referencing representative entries."""
        metadata_path = tmp_path / "metadata.tsv"
        # Existing file has GCF_A02 as SpA representative (not first alpha)
        self._write_metadata(metadata_path, include_representative=True)
        output_path = tmp_path / "output.tsv"

        result = runner.invoke(
            app,
            [
                "select-representatives",
                "--metadata", str(metadata_path),
                "--output", str(output_path),
            ],
        )

        assert result.exit_code == 0

        df = pl.read_csv(output_path, separator="\t")
        # GCF_A02 was the existing self-referencing rep for SpA
        spa_reps = df.filter(pl.col("species") == "SpA")["representative"].to_list()
        assert all(r == "GCF_A02" for r in spa_reps)

    def test_success_message(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should print summary with species and genome counts."""
        metadata_path = tmp_path / "metadata.tsv"
        self._write_metadata(metadata_path)
        output_path = tmp_path / "output.tsv"

        result = runner.invoke(
            app,
            [
                "select-representatives",
                "--metadata", str(metadata_path),
                "--output", str(output_path),
            ],
        )

        assert result.exit_code == 0
        assert "2" in result.output  # 2 species or representatives
        assert "4" in result.output  # 4 genomes

    def test_missing_required_columns(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should exit with error if required columns are missing."""
        bad_metadata = tmp_path / "bad_metadata.tsv"
        pl.DataFrame({"accession": ["A"], "genus": ["G"]}).write_csv(
            bad_metadata, separator="\t"
        )
        output_path = tmp_path / "output.tsv"

        result = runner.invoke(
            app,
            [
                "select-representatives",
                "--metadata", str(bad_metadata),
                "--output", str(output_path),
            ],
        )

        assert result.exit_code == 1

    def test_quiet_mode(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should suppress output in quiet mode."""
        metadata_path = tmp_path / "metadata.tsv"
        self._write_metadata(metadata_path)
        output_path = tmp_path / "output.tsv"

        result = runner.invoke(
            app,
            [
                "select-representatives",
                "--metadata", str(metadata_path),
                "--output", str(output_path),
                "--quiet",
            ],
        )

        assert result.exit_code == 0

    def test_creates_output_directory(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should create parent directories for output path."""
        metadata_path = tmp_path / "metadata.tsv"
        self._write_metadata(metadata_path)
        output_path = tmp_path / "nested" / "deep" / "output.tsv"

        result = runner.invoke(
            app,
            [
                "select-representatives",
                "--metadata", str(metadata_path),
                "--output", str(output_path),
            ],
        )

        assert result.exit_code == 0
        assert output_path.exists()

    def test_representative_column_position(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Representative column should appear after family column."""
        metadata_path = tmp_path / "metadata.tsv"
        self._write_metadata(metadata_path)
        output_path = tmp_path / "output.tsv"

        result = runner.invoke(
            app,
            [
                "select-representatives",
                "--metadata", str(metadata_path),
                "--output", str(output_path),
            ],
        )

        assert result.exit_code == 0
        df = pl.read_csv(output_path, separator="\t")
        cols = df.columns
        family_idx = cols.index("family")
        rep_idx = cols.index("representative")
        assert rep_idx == family_idx + 1
