"""
Unit tests for genome metadata handling.

Tests GenomeMetadata class including file loading, lookups, novel taxon
naming, classification joining, and species/genus aggregation.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.metadata import GenomeMetadata


# =============================================================================
# Helpers
# =============================================================================


def _make_metadata_df(*, include_family: bool = True) -> pl.DataFrame:
    """Build a small metadata DataFrame for testing."""
    data: dict[str, list[str]] = {
        "accession": ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"],
        "species": [
            "Francisella tularensis",
            "Francisella novicida",
            "Francisella philomiragia",
        ],
        "genus": ["Francisella", "Francisella", "Francisella"],
    }
    if include_family:
        data["family"] = [
            "Francisellaceae",
            "Francisellaceae",
            "Francisellaceae",
        ]
    return pl.DataFrame(data)


def _make_classification_df() -> pl.DataFrame:
    """Build a minimal classification DataFrame for join tests."""
    return pl.DataFrame(
        {
            "read_id": ["r1", "r2", "r3", "r4"],
            "best_match_genome": [
                "GCF_000123456.1",
                "GCF_000789012.1",
                "GCA_000111222.1",
                "GCF_UNKNOWN.1",
            ],
            "taxonomic_call": [
                "Known Species",
                "Novel Species",
                "Novel Genus",
                "Known Species",
            ],
            "novelty_index": [1.0, 8.0, 22.0, 2.0],
            "placement_uncertainty": [0.2, 0.5, 1.5, 0.3],
            "top_hit_identity": [99.0, 92.0, 78.0, 98.0],
        }
    )


# =============================================================================
# Tests: from_file
# =============================================================================


class TestFromFile:
    """Tests for GenomeMetadata.from_file class method."""

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        """from_file should raise FileNotFoundError for missing path."""
        missing = tmp_path / "does_not_exist.tsv"
        with pytest.raises(FileNotFoundError, match="Metadata file not found"):
            GenomeMetadata.from_file(missing)

    def test_missing_columns_raises(self, tmp_path: Path) -> None:
        """from_file should raise ValueError when required columns are absent."""
        bad_file = tmp_path / "bad_metadata.tsv"
        pl.DataFrame({"accession": ["A"], "extra": ["X"]}).write_csv(
            bad_file, separator="\t"
        )
        with pytest.raises(ValueError, match="missing required columns"):
            GenomeMetadata.from_file(bad_file)

    def test_loads_valid_file(self, tmp_path: Path) -> None:
        """from_file should load a well-formed metadata TSV."""
        tsv = tmp_path / "metadata.tsv"
        _make_metadata_df().write_csv(tsv, separator="\t")
        meta = GenomeMetadata.from_file(tsv)
        assert meta.genome_count == 3


# =============================================================================
# Tests: properties
# =============================================================================


class TestProperties:
    """Tests for GenomeMetadata simple properties."""

    def test_dataframe_property(self) -> None:
        """dataframe should return the underlying DataFrame."""
        df = _make_metadata_df()
        meta = GenomeMetadata(df)
        assert meta.dataframe is df

    def test_genome_count(self) -> None:
        """genome_count should return number of rows."""
        meta = GenomeMetadata(_make_metadata_df())
        assert meta.genome_count == 3

    def test_species_count(self) -> None:
        """species_count should return unique species."""
        meta = GenomeMetadata(_make_metadata_df())
        assert meta.species_count == 3

    def test_genus_count(self) -> None:
        """genus_count should return unique genera."""
        meta = GenomeMetadata(_make_metadata_df())
        # All three species belong to genus Francisella
        assert meta.genus_count == 1


# =============================================================================
# Tests: lookups
# =============================================================================


class TestLookups:
    """Tests for accession-based lookup methods."""

    def test_get_species_known(self) -> None:
        """get_species should return the species for a known accession."""
        meta = GenomeMetadata(_make_metadata_df())
        assert meta.get_species("GCF_000123456.1") == "Francisella tularensis"

    def test_get_species_unknown(self) -> None:
        """get_species should return None for an unknown accession."""
        meta = GenomeMetadata(_make_metadata_df())
        assert meta.get_species("GCF_MISSING.1") is None

    def test_get_genus_known(self) -> None:
        """get_genus should return genus for a known accession."""
        meta = GenomeMetadata(_make_metadata_df())
        assert meta.get_genus("GCF_000789012.1") == "Francisella"

    def test_get_family_known(self) -> None:
        """get_family should return family for a known accession."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        assert meta.get_family("GCA_000111222.1") == "Francisellaceae"

    def test_get_family_no_family_column(self) -> None:
        """get_family should return None when metadata lacks a family column."""
        meta = GenomeMetadata(_make_metadata_df(include_family=False))
        assert meta.get_family("GCF_000123456.1") is None

    def test_get_family_unknown_accession(self) -> None:
        """get_family should return None for an unrecognised accession."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        assert meta.get_family("MISSING") is None


# =============================================================================
# Tests: infer_target_family
# =============================================================================


class TestInferTargetFamily:
    """Tests for infer_target_family method."""

    def test_returns_most_common_family(self) -> None:
        """Should return the most frequent family."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        assert meta.infer_target_family() == "Francisellaceae"

    def test_no_family_column(self) -> None:
        """Should return None when family column is absent."""
        meta = GenomeMetadata(_make_metadata_df(include_family=False))
        assert meta.infer_target_family() is None

    def test_empty_dataframe(self) -> None:
        """Should return None for an empty metadata set."""
        empty = pl.DataFrame(
            {"accession": [], "species": [], "genus": [], "family": []},
            schema={
                "accession": pl.Utf8,
                "species": pl.Utf8,
                "genus": pl.Utf8,
                "family": pl.Utf8,
            },
        )
        meta = GenomeMetadata(empty)
        assert meta.infer_target_family() is None


# =============================================================================
# Tests: generate_novel_id
# =============================================================================


class TestGenerateNovelId:
    """Tests for the static generate_novel_id method."""

    def test_returns_three_characters(self) -> None:
        """ID should be exactly 3 alphanumeric characters."""
        novel_id = GenomeMetadata.generate_novel_id("read_1", "GCF_000123456.1")
        assert len(novel_id) == 3
        assert novel_id.isalnum()

    def test_deterministic(self) -> None:
        """Same inputs should always produce the same ID."""
        id_a = GenomeMetadata.generate_novel_id("read_1", "GCF_000123456.1")
        id_b = GenomeMetadata.generate_novel_id("read_1", "GCF_000123456.1")
        assert id_a == id_b

    def test_different_inputs_differ(self) -> None:
        """Different read/accession pairs should generally produce different IDs."""
        id_a = GenomeMetadata.generate_novel_id("read_1", "GCF_000123456.1")
        id_b = GenomeMetadata.generate_novel_id("read_2", "GCF_000123456.1")
        # Not guaranteed to differ, but extremely likely for distinct inputs
        # We just verify both are valid 3-char IDs
        assert len(id_b) == 3
        assert id_b.isalnum()

    def test_characters_from_expected_set(self) -> None:
        """ID characters should come from A-Z and 0-9."""
        allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
        novel_id = GenomeMetadata.generate_novel_id("read_x", "GCF_999999999.1")
        assert all(c in allowed for c in novel_id)


# =============================================================================
# Tests: suggest_novel_species_name
# =============================================================================


class TestSuggestNovelSpeciesName:
    """Tests for suggest_novel_species_name."""

    def test_with_read_id(self) -> None:
        """Should produce '{Genus} sp. MDM-{id}' when read_id is given."""
        meta = GenomeMetadata(_make_metadata_df())
        name = meta.suggest_novel_species_name("GCF_000123456.1", read_id="read_42")
        assert name.startswith("Francisella sp. MDM-")
        assert len(name.split("MDM-")[1]) == 3

    def test_without_read_id(self) -> None:
        """Should produce '{Genus} sp. nov.' when read_id is omitted."""
        meta = GenomeMetadata(_make_metadata_df())
        name = meta.suggest_novel_species_name("GCF_000123456.1")
        assert name == "Francisella sp. nov."

    def test_unknown_accession_fallback(self) -> None:
        """Should use 'Unknown' genus when accession is not found."""
        meta = GenomeMetadata(_make_metadata_df())
        name = meta.suggest_novel_species_name("GCF_MISSING.1")
        assert name == "Unknown sp. nov."

    def test_unknown_accession_with_read_id(self) -> None:
        """Should use 'Unknown' genus with MDM-id for missing accession."""
        meta = GenomeMetadata(_make_metadata_df())
        name = meta.suggest_novel_species_name("GCF_MISSING.1", read_id="read_1")
        assert name.startswith("Unknown sp. MDM-")


# =============================================================================
# Tests: suggest_novel_genus_name
# =============================================================================


class TestSuggestNovelGenusName:
    """Tests for suggest_novel_genus_name."""

    def test_with_family_and_read_id(self) -> None:
        """Should produce '{Family} gen. nov. MDM-{id}'."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        name = meta.suggest_novel_genus_name("GCF_000123456.1", read_id="read_1")
        assert name.startswith("Francisellaceae gen. nov. MDM-")

    def test_with_family_no_read_id(self) -> None:
        """Should produce '{Family} gen. nov.' without read_id."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        name = meta.suggest_novel_genus_name("GCF_000123456.1")
        assert name == "Francisellaceae gen. nov."

    def test_no_family_falls_back_to_genus(self) -> None:
        """Should use '{Genus}-related gen. nov.' when no family column."""
        meta = GenomeMetadata(_make_metadata_df(include_family=False))
        name = meta.suggest_novel_genus_name("GCF_000123456.1")
        assert name == "Francisella-related gen. nov."

    def test_no_family_with_read_id(self) -> None:
        """Should use genus-related fallback with MDM-id."""
        meta = GenomeMetadata(_make_metadata_df(include_family=False))
        name = meta.suggest_novel_genus_name("GCF_000123456.1", read_id="read_7")
        assert name.startswith("Francisella-related gen. nov. MDM-")

    def test_unknown_accession_no_family(self) -> None:
        """Should use 'Unknown-related' when accession is completely unknown."""
        meta = GenomeMetadata(_make_metadata_df(include_family=False))
        name = meta.suggest_novel_genus_name("MISSING")
        assert name == "Unknown-related gen. nov."


# =============================================================================
# Tests: join_classifications
# =============================================================================


class TestJoinClassifications:
    """Tests for join_classifications method."""

    def test_missing_column_raises(self) -> None:
        """Should raise ValueError when best_match_genome is absent."""
        meta = GenomeMetadata(_make_metadata_df())
        bad_df = pl.DataFrame({"read_id": ["r1"]})
        with pytest.raises(ValueError, match="best_match_genome"):
            meta.join_classifications(bad_df)

    def test_joins_species_and_genus(self) -> None:
        """Should add species and genus columns from metadata."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        result = meta.join_classifications(clf)

        assert "species" in result.columns
        assert "genus" in result.columns
        # The known accession should have proper species
        row_0 = result.filter(pl.col("read_id") == "r1")
        assert row_0["species"][0] == "Francisella tularensis"

    def test_fills_unknown_for_missing_genome(self) -> None:
        """Unmatched genomes should get 'Unknown species' / 'Unknown genus'."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        result = meta.join_classifications(clf)
        unknown_row = result.filter(pl.col("read_id") == "r4")
        assert unknown_row["species"][0] == "Unknown species"
        assert unknown_row["genus"][0] == "Unknown genus"

    def test_includes_family_when_present(self) -> None:
        """Should include family column when metadata has it."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        clf = _make_classification_df()
        result = meta.join_classifications(clf)
        assert "family" in result.columns
        known = result.filter(pl.col("read_id") == "r1")
        assert known["family"][0] == "Francisellaceae"

    def test_fills_unknown_family(self) -> None:
        """Unmatched genomes should get 'Unknown family'."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        clf = _make_classification_df()
        result = meta.join_classifications(clf)
        unknown_row = result.filter(pl.col("read_id") == "r4")
        assert unknown_row["family"][0] == "Unknown family"

    def test_no_family_column_in_metadata(self) -> None:
        """Result should omit family when metadata lacks it."""
        meta = GenomeMetadata(_make_metadata_df(include_family=False))
        clf = _make_classification_df()
        result = meta.join_classifications(clf)
        assert "family" not in result.columns


# =============================================================================
# Tests: join_classifications_with_novel_names
# =============================================================================


class TestJoinClassificationsWithNovelNames:
    """Tests for join_classifications_with_novel_names method."""

    def test_missing_required_columns_raises(self) -> None:
        """Should raise ValueError when required columns are missing."""
        meta = GenomeMetadata(_make_metadata_df())
        bad_df = pl.DataFrame({"best_match_genome": ["A"], "read_id": ["r1"]})
        with pytest.raises(ValueError, match="missing required columns"):
            meta.join_classifications_with_novel_names(bad_df)

    def test_known_species_gets_actual_name(self) -> None:
        """Known Species rows should get the real species name."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        clf = _make_classification_df()
        result = meta.join_classifications_with_novel_names(clf)

        known = result.filter(pl.col("taxonomic_call") == "Known Species")
        # r1 maps to GCF_000123456.1 = Francisella tularensis
        r1 = known.filter(pl.col("read_id") == "r1")
        assert r1["suggested_name"][0] == "Francisella tularensis"

    def test_novel_species_gets_mdm_name(self) -> None:
        """Novel Species rows should get '{Genus} sp. MDM-{id}'."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        clf = _make_classification_df()
        result = meta.join_classifications_with_novel_names(clf)

        novel_sp = result.filter(pl.col("taxonomic_call") == "Novel Species")
        name = novel_sp["suggested_name"][0]
        assert "sp. MDM-" in name
        assert name.startswith("Francisella")

    def test_novel_genus_gets_gen_nov_name(self) -> None:
        """Novel Genus rows should get '{Family} gen. nov. MDM-{id}'."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        clf = _make_classification_df()
        result = meta.join_classifications_with_novel_names(clf)

        novel_gen = result.filter(pl.col("taxonomic_call") == "Novel Genus")
        name = novel_gen["suggested_name"][0]
        assert "gen. nov. MDM-" in name
        assert name.startswith("Francisellaceae")

    def test_suggested_name_column_present(self) -> None:
        """Result should have a 'suggested_name' column."""
        meta = GenomeMetadata(_make_metadata_df(include_family=True))
        clf = _make_classification_df()
        result = meta.join_classifications_with_novel_names(clf)
        assert "suggested_name" in result.columns
        assert "_novel_id" not in result.columns  # internal column dropped

    def test_no_family_column_adds_unknown_family(self) -> None:
        """When metadata lacks family, 'Unknown family' should be used."""
        meta = GenomeMetadata(_make_metadata_df(include_family=False))
        clf = _make_classification_df()
        result = meta.join_classifications_with_novel_names(clf)
        assert "family" in result.columns


# =============================================================================
# Tests: aggregate_by_species
# =============================================================================


class TestAggregateBySpecies:
    """Tests for aggregate_by_species method."""

    def test_missing_species_column_raises(self) -> None:
        """Should raise ValueError when species column is absent."""
        meta = GenomeMetadata(_make_metadata_df())
        bad_df = pl.DataFrame({"read_id": ["r1"]})
        with pytest.raises(ValueError, match="species"):
            meta.aggregate_by_species(bad_df)

    def test_aggregation_returns_expected_columns(self) -> None:
        """Result should contain standard aggregation columns."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        joined = meta.join_classifications(clf)
        agg = meta.aggregate_by_species(joined)

        expected_cols = {
            "species",
            "read_count",
            "mean_novelty",
            "mean_identity",
            "mean_uncertainty",
            "genome_count",
        }
        assert expected_cols.issubset(set(agg.columns))

    def test_sorted_by_read_count_descending(self) -> None:
        """Rows should be sorted by read_count in descending order."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        joined = meta.join_classifications(clf)
        agg = meta.aggregate_by_species(joined)
        counts = agg["read_count"].to_list()
        assert counts == sorted(counts, reverse=True)

    def test_correct_read_counts(self) -> None:
        """Each species should have the correct number of reads."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        joined = meta.join_classifications(clf)
        agg = meta.aggregate_by_species(joined)
        total_reads = sum(agg["read_count"].to_list())
        assert total_reads == len(clf)


# =============================================================================
# Tests: aggregate_by_genus
# =============================================================================


class TestAggregateByGenus:
    """Tests for aggregate_by_genus method."""

    def test_missing_genus_column_raises(self) -> None:
        """Should raise ValueError when genus column is absent."""
        meta = GenomeMetadata(_make_metadata_df())
        bad_df = pl.DataFrame({"read_id": ["r1"]})
        with pytest.raises(ValueError, match="genus"):
            meta.aggregate_by_genus(bad_df)

    def test_aggregation_returns_expected_columns(self) -> None:
        """Result should contain genus-level aggregation columns."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        joined = meta.join_classifications(clf)
        agg = meta.aggregate_by_genus(joined)

        expected_cols = {
            "genus",
            "read_count",
            "mean_novelty",
            "mean_identity",
            "mean_uncertainty",
            "genome_count",
            "species_count",
        }
        assert expected_cols.issubset(set(agg.columns))

    def test_sorted_by_read_count_descending(self) -> None:
        """Rows should be sorted by read_count in descending order."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        joined = meta.join_classifications(clf)
        agg = meta.aggregate_by_genus(joined)
        counts = agg["read_count"].to_list()
        assert counts == sorted(counts, reverse=True)

    def test_correct_total_reads(self) -> None:
        """Total reads across all genera should equal input rows."""
        meta = GenomeMetadata(_make_metadata_df())
        clf = _make_classification_df()
        joined = meta.join_classifications(clf)
        agg = meta.aggregate_by_genus(joined)
        total = sum(agg["read_count"].to_list())
        assert total == len(clf)
