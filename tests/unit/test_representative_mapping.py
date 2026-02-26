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


# =============================================================================
# 8. Family validation with representative_mapping (broad database scenario)
# =============================================================================


class TestFamilyValidationWithRepresentativeMapping:
    """Tests for family validation when using representative_mapping.

    Scenario: broad alignment database with GTDB representatives for all
    bacteria plus all genomes within the focus family. The ANI matrix
    contains only family representatives, but the representative_mapping
    keys cover all family genomes. Non-family genomes (from the broad
    database) should be treated as external, while non-representative
    family genomes should be treated as in-family.
    """

    @pytest.fixture
    def broad_db_setup(self) -> dict:
        """Set up a broad-database scenario with family and non-family genomes.

        Family genomes: GCF_F01 (rep), GCF_F02, GCF_F03 (rep), GCF_F04
        Non-family genomes: GCF_EXT1, GCF_EXT2 (not in mapping)
        ANI matrix: GCF_F01 and GCF_F03 only (representatives)
        """
        representative_mapping = {
            "GCF_F01": "GCF_F01",  # rep for species A
            "GCF_F02": "GCF_F01",  # non-rep, maps to GCF_F01
            "GCF_F03": "GCF_F03",  # rep for species B
            "GCF_F04": "GCF_F03",  # non-rep, maps to GCF_F03
        }

        ani_dict = {
            "GCF_F01": {"GCF_F03": 82.0},
            "GCF_F03": {"GCF_F01": 82.0},
        }

        # BLAST hits: reads hit both family and non-family genomes.
        # Family genomes get higher bitscores so reads are not off-target.
        blast_rows = []
        for read_idx in range(3):
            # Hits to family genomes (high bitscore)
            for genome, bitscore, pident in [
                ("GCF_F01", 300.0, 98.0),
                ("GCF_F02", 295.0, 97.5),
                ("GCF_F03", 280.0, 95.0),
                ("GCF_F04", 275.0, 94.5),
            ]:
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}|contig_1",
                    "pident": pident,
                    "length": 200,
                    "mismatch": int(200 * (1 - pident / 100)),
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 200,
                    "sstart": 500,
                    "send": 700,
                    "evalue": 1e-40,
                    "bitscore": bitscore,
                    "genome_name": genome,
                })
            # Hits to non-family genomes (lower bitscore)
            for genome, bitscore, pident in [
                ("GCF_EXT1", 200.0, 85.0),
                ("GCF_EXT2", 180.0, 82.0),
            ]:
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}|contig_1",
                    "pident": pident,
                    "length": 200,
                    "mismatch": int(200 * (1 - pident / 100)),
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 200,
                    "sstart": 500,
                    "send": 700,
                    "evalue": 1e-20,
                    "bitscore": bitscore,
                    "genome_name": genome,
                })

        blast_df = pl.DataFrame(blast_rows)

        return {
            "representative_mapping": representative_mapping,
            "ani_dict": ani_dict,
            "blast_df": blast_df,
        }

    def test_non_rep_family_genomes_treated_as_in_family(
        self, broad_db_setup: dict
    ) -> None:
        """Non-representative family genomes should be in-family, not external.

        When representative_mapping is provided, its keys define the set of
        in-family genomes. GCF_F02 and GCF_F04 are in the mapping so they
        should be treated as in-family even though they are not in the ANI
        matrix directly.
        """
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )
        from metadarkmatter.models.config import ScoringConfig

        ani_matrix = ANIMatrix(broad_db_setup["ani_dict"])
        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            config=config,
            representative_mapping=broad_db_setup["representative_mapping"],
        )

        result = classifier.classify_file(broad_db_setup["blast_df"])

        # No reads should be Off-target because all reads have stronger
        # in-family hits than external hits
        off_target = result.filter(pl.col("taxonomic_call") == "Off-target")
        assert len(off_target) == 0, (
            f"Expected no Off-target reads, got {len(off_target)}. "
            "Non-rep family genomes may be incorrectly treated as external."
        )

    def test_non_family_genomes_treated_as_external(
        self, broad_db_setup: dict
    ) -> None:
        """Non-family genomes (not in representative_mapping) should be external.

        GCF_EXT1 and GCF_EXT2 are not keys in the representative_mapping,
        so they should be treated as external hits for family validation.
        """
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )
        from metadarkmatter.models.config import ScoringConfig

        ani_matrix = ANIMatrix(broad_db_setup["ani_dict"])
        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            config=config,
            representative_mapping=broad_db_setup["representative_mapping"],
        )

        result = classifier.classify_file(broad_db_setup["blast_df"])

        # Family validation columns should exist (family validation was active)
        assert "family_bitscore_ratio" in result.columns
        assert "in_family_hit_fraction" in result.columns

        # in_family_hit_fraction should reflect 4 family hits out of 6 total
        fractions = result["in_family_hit_fraction"].to_list()
        expected_fraction = 4.0 / 6.0
        for frac in fractions:
            assert abs(frac - expected_fraction) < 0.01, (
                f"Expected in_family_hit_fraction ~{expected_fraction:.3f}, got {frac}"
            )

    def test_fallback_to_ani_matrix_without_mapping(
        self, broad_db_setup: dict
    ) -> None:
        """Without representative_mapping, in-family is defined by ANI matrix membership.

        This tests backwards compatibility: when no representative_mapping is
        provided, the original behavior (ANI matrix genomes = in-family) applies.
        """
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )
        from metadarkmatter.models.config import ScoringConfig

        ani_matrix = ANIMatrix(broad_db_setup["ani_dict"])
        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            config=config,
            representative_mapping=None,  # No mapping
        )

        result = classifier.classify_file(broad_db_setup["blast_df"])

        # Without mapping, only GCF_F01 and GCF_F03 (ANI matrix genomes)
        # are in-family. GCF_F02, GCF_F04, GCF_EXT1, GCF_EXT2 are all external.
        # in_family_hit_fraction should be 2/6
        fractions = result["in_family_hit_fraction"].to_list()
        expected_fraction = 2.0 / 6.0
        for frac in fractions:
            assert abs(frac - expected_fraction) < 0.01, (
                f"Without mapping, expected in_family_hit_fraction ~{expected_fraction:.3f}, got {frac}"
            )

    def test_off_target_when_external_hits_dominate(self) -> None:
        """Reads with better external hits should be classified as Off-target.

        Create a scenario where external genomes have higher bitscores
        than family genomes, triggering the off-target threshold.
        """
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )
        from metadarkmatter.models.config import ScoringConfig

        representative_mapping = {
            "GCF_F01": "GCF_F01",
        }
        ani_dict = {
            "GCF_F01": {"GCF_F01": 100.0},
        }

        # Read has a weak family hit and a strong external hit
        blast_df = pl.DataFrame([
            {
                "qseqid": "read_001",
                "sseqid": "GCF_F01|contig_1",
                "pident": 80.0,
                "length": 200,
                "mismatch": 40,
                "gapopen": 0,
                "qstart": 1,
                "qend": 200,
                "sstart": 500,
                "send": 700,
                "evalue": 1e-20,
                "bitscore": 100.0,
                "genome_name": "GCF_F01",
            },
            {
                "qseqid": "read_001",
                "sseqid": "GCF_EXT1|contig_1",
                "pident": 99.0,
                "length": 200,
                "mismatch": 2,
                "gapopen": 0,
                "qstart": 1,
                "qend": 200,
                "sstart": 500,
                "send": 700,
                "evalue": 1e-80,
                "bitscore": 350.0,
                "genome_name": "GCF_EXT1",
            },
        ])

        ani_matrix = ANIMatrix(ani_dict)
        config = ScoringConfig(target_family="f__TestFamily", family_ratio_threshold=0.8)
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            config=config,
            representative_mapping=representative_mapping,
        )

        result = classifier.classify_file(blast_df)

        # bitscore ratio = 100/350 = 0.286, well below 0.8 threshold
        assert result["taxonomic_call"][0] == "Off-target"


# =============================================================================
# 9. validate_ani_genome_coverage with broad database
# =============================================================================


class TestValidateAniGenomeCoverageBroadDatabase:
    """Tests for validate_ani_genome_coverage with broad databases.

    When using a broad alignment database (all bacteria + focus family),
    the BLAST file contains hits to both family and non-family genomes.
    Only family genomes (those in representative_mapping) should be
    checked against the ANI matrix.
    """

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

    def test_broad_db_excludes_non_family_genomes(self, tmp_path: Path) -> None:
        """Non-family genomes should be excluded from coverage calculation.

        When representative_mapping is provided, only genomes present as
        keys in the mapping (family genomes) should be checked. External
        genomes from the broad database should not inflate the denominator.
        """
        from metadarkmatter.cli.score import validate_ani_genome_coverage

        ani_dict = {
            "GCF_REP1": {"GCF_REP2": 90.0},
            "GCF_REP2": {"GCF_REP1": 90.0},
        }
        ani_matrix = ANIMatrix(ani_dict)

        # BLAST file: 2 family genomes + 3 non-family genomes
        blast_path = tmp_path / "blast.tsv"
        self._write_blast_file(
            blast_path,
            ["GCF_FAM1", "GCF_FAM2", "GCF_EXT1", "GCF_EXT2", "GCF_EXT3"],
        )

        # Only family genomes are in the mapping
        rep_mapping = {
            "GCF_FAM1": "GCF_REP1",
            "GCF_FAM2": "GCF_REP2",
        }

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix, representative_mapping=rep_mapping,
        )

        # Only 2 family genomes checked, both map to representatives in ANI
        assert total == 2, f"Expected 2 family genomes, got {total}"
        assert matched == 2
        assert pct == 100.0
        assert len(missing) == 0

    def test_broad_db_without_mapping_includes_all(self, tmp_path: Path) -> None:
        """Without mapping, all genomes are checked (backwards compatible).

        This verifies that the fix does not break the no-mapping case.
        """
        from metadarkmatter.cli.score import validate_ani_genome_coverage

        ani_dict = {
            "GCF_REP1": {"GCF_REP2": 90.0},
            "GCF_REP2": {"GCF_REP1": 90.0},
        }
        ani_matrix = ANIMatrix(ani_dict)

        blast_path = tmp_path / "blast.tsv"
        self._write_blast_file(
            blast_path,
            ["GCF_FAM1", "GCF_FAM2", "GCF_EXT1"],
        )

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix, representative_mapping=None,
        )

        # All 3 genomes checked, none in ANI matrix
        assert total == 3
        assert matched == 0
        assert pct == 0.0

    def test_broad_db_empty_blast_genomes_in_mapping(self, tmp_path: Path) -> None:
        """When no BLAST genomes are in the mapping, coverage should be 0/0."""
        from metadarkmatter.cli.score import validate_ani_genome_coverage

        ani_dict = {
            "GCF_REP1": {"GCF_REP1": 100.0},
        }
        ani_matrix = ANIMatrix(ani_dict)

        # BLAST file has only non-family genomes
        blast_path = tmp_path / "blast.tsv"
        self._write_blast_file(blast_path, ["GCF_EXT1", "GCF_EXT2"])

        rep_mapping = {
            "GCF_FAM1": "GCF_REP1",  # No BLAST hits to family genomes
        }

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix, representative_mapping=rep_mapping,
        )

        # No family genomes in BLAST file
        assert total == 0
        assert matched == 0
        assert pct == 0.0


# =============================================================================
# 10. Representative collapse bug fix
# =============================================================================


class TestRepresentativeCollapsePlacementUncertainty:
    """Tests for placement_uncertainty when all hits collapse to one representative.

    Bug scenario: a read hits 10 strains of the same species, all mapping to
    the same representative. Before the fix, num_secondary_genomes == 0 but
    num_ambiguous_hits > 1, causing a left-join to produce ANI = 0 and
    placement_uncertainty = 100%. After the fix, reads with
    num_secondary_genomes == 0 should get placement_uncertainty = 0.
    """

    @pytest.fixture
    def single_rep_setup(self) -> dict:
        """Set up a scenario where all hits collapse to one representative.

        4 genomes (strains), all mapping to one representative (GCF_REP1).
        ANI matrix has only GCF_REP1. Reads hit all 4 strains with high identity.
        """
        representative_mapping = {
            "GCF_S01": "GCF_REP1",
            "GCF_S02": "GCF_REP1",
            "GCF_S03": "GCF_REP1",
            "GCF_S04": "GCF_REP1",
        }

        # ANI matrix with only the representative (single-genome matrix)
        ani_dict = {
            "GCF_REP1": {},
        }

        # BLAST hits: reads hit all 4 strains with very high identity
        blast_rows = []
        for read_idx in range(3):
            for strain_idx, genome in enumerate(
                ["GCF_S01", "GCF_S02", "GCF_S03", "GCF_S04"]
            ):
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}|contig_1",
                    "pident": 99.3 - strain_idx * 0.1,
                    "length": 250,
                    "mismatch": 2 + strain_idx,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 250,
                    "sstart": 1000,
                    "send": 1250,
                    "evalue": 1e-80,
                    "bitscore": 480.0 - strain_idx * 2.0,
                    "genome_name": genome,
                })

        blast_df = pl.DataFrame(blast_rows)

        return {
            "representative_mapping": representative_mapping,
            "ani_dict": ani_dict,
            "blast_df": blast_df,
        }

    def test_collapsed_reads_get_zero_uncertainty(self, single_rep_setup: dict) -> None:
        """Reads where all hits collapse to one representative should get 0% uncertainty."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(single_rep_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=single_rep_setup["representative_mapping"],
        )

        result = classifier.classify_file(single_rep_setup["blast_df"])

        assert len(result) == 3
        for u in result["placement_uncertainty"].to_list():
            assert u == 0.0, (
                f"Expected 0% uncertainty for collapsed reads, got {u}%"
            )

    def test_collapsed_reads_classified_as_known(self, single_rep_setup: dict) -> None:
        """High-identity reads with collapsed representatives should be Known Species."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(single_rep_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=single_rep_setup["representative_mapping"],
        )

        result = classifier.classify_file(single_rep_setup["blast_df"])

        for call in result["taxonomic_call"].to_list():
            assert call == "Known Species", (
                f"Expected Known Species for high-identity collapsed reads, got {call}"
            )

    def test_collapsed_reads_get_inferred_uncertainty(
        self, single_rep_setup: dict
    ) -> None:
        """Collapsed reads should get inferred uncertainty (like single-hit reads)."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(single_rep_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=single_rep_setup["representative_mapping"],
        )

        result = classifier.classify_file(single_rep_setup["blast_df"])

        for utype in result["uncertainty_type"].to_list():
            assert utype == "inferred", (
                f"Expected inferred uncertainty for collapsed reads, got {utype}"
            )
        # Inferred uncertainty should be non-null
        for iu in result["inferred_uncertainty"].to_list():
            assert iu is not None, "inferred_uncertainty should not be null"

    def test_collapsed_reads_ambiguity_scope_unambiguous(
        self, single_rep_setup: dict
    ) -> None:
        """Collapsed reads should have ambiguity_scope = 'unambiguous'."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        ani_matrix = ANIMatrix(single_rep_setup["ani_dict"])
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=single_rep_setup["representative_mapping"],
        )

        result = classifier.classify_file(single_rep_setup["blast_df"])

        for scope in result["ambiguity_scope"].to_list():
            assert scope == "unambiguous", (
                f"Expected unambiguous scope for collapsed reads, got {scope}"
            )


class TestMultiRepresentativeStillWorks:
    """Regression test: reads hitting multiple distinct representatives still get ANI-based uncertainty."""

    def test_multi_rep_reads_use_ani_uncertainty(self) -> None:
        """Reads hitting genomes from different species should still use ANI for uncertainty."""
        from metadarkmatter.core.classification.classifiers.vectorized import (
            VectorizedClassifier,
        )

        # Two species, each with strains
        representative_mapping = {
            "GCF_A1": "GCF_REP_A",
            "GCF_A2": "GCF_REP_A",
            "GCF_B1": "GCF_REP_B",
            "GCF_B2": "GCF_REP_B",
        }

        ani_dict = {
            "GCF_REP_A": {"GCF_REP_B": 85.0},
            "GCF_REP_B": {"GCF_REP_A": 85.0},
        }

        blast_rows = []
        for genome, pident, bitscore in [
            ("GCF_A1", 98.0, 500),
            ("GCF_A2", 97.5, 495),
            ("GCF_B1", 96.0, 490),
            ("GCF_B2", 95.5, 485),
        ]:
            blast_rows.append({
                "qseqid": "read_001",
                "sseqid": f"{genome}|contig_1",
                "pident": pident,
                "length": 250,
                "mismatch": int(250 * (1 - pident / 100)),
                "gapopen": 0,
                "qstart": 1,
                "qend": 250,
                "sstart": 1000,
                "send": 1250,
                "evalue": 1e-80,
                "bitscore": float(bitscore),
                "genome_name": genome,
            })

        blast_df = pl.DataFrame(blast_rows)

        ani_matrix = ANIMatrix(ani_dict)
        classifier = VectorizedClassifier(
            ani_matrix=ani_matrix,
            representative_mapping=representative_mapping,
        )

        result = classifier.classify_file(blast_df)

        assert len(result) == 1
        u = result["placement_uncertainty"][0]
        # ANI between reps is 85%, so uncertainty = 100 - 85 = 15%
        assert u == pytest.approx(15.0, abs=1.0), (
            f"Expected ~15% uncertainty for multi-rep reads, got {u}%"
        )
        assert result["uncertainty_type"][0] == "measured"
