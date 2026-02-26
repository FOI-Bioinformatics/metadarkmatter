"""
Unit tests for novel diversity clustering.

Tests NovelCluster and NovelDiversitySummary models, as well as the
NovelDiversityAnalyzer clustering logic.
"""

from __future__ import annotations

import pytest
import polars as pl

from metadarkmatter.core.novel_diversity import (
    NovelCluster,
    NovelDiversityAnalyzer,
    NovelDiversitySummary,
)


# =============================================================================
# NovelCluster Model Tests
# =============================================================================


class TestNovelCluster:
    """Tests for NovelCluster Pydantic model."""

    @pytest.fixture
    def novel_species_cluster(self):
        """Sample Novel Species cluster."""
        return NovelCluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_000123456.1",
            nearest_species="Francisella tularensis",
            nearest_genus="Francisella",
            nearest_family="Francisellaceae",
            novelty_band=5,
            read_count=25,
            mean_novelty_index=7.5,
            novelty_min=5.2,
            novelty_max=9.8,
            mean_placement_uncertainty=1.2,
            mean_discovery_score=82.5,
            suggested_name="Francisella sp. MDM-001",
            confidence="High",
            phylogenetic_context="Novel species within Francisella",
            read_ids=["read_001", "read_002", "read_003"],
        )

    @pytest.fixture
    def novel_genus_cluster(self):
        """Sample Novel Genus cluster."""
        return NovelCluster(
            cluster_id="NGN_001",
            taxonomic_call="Novel Genus",
            nearest_genome="GCF_000789012.1",
            nearest_species="Francisella novicida",
            nearest_genus="Francisella",
            nearest_family="Francisellaceae",
            novelty_band=20,
            read_count=12,
            mean_novelty_index=22.3,
            novelty_min=20.1,
            novelty_max=24.5,
            mean_placement_uncertainty=3.5,
            mean_discovery_score=65.0,
            suggested_name="Francisellaceae gen. nov. MDM-001",
            confidence="Medium",
            phylogenetic_context="Novel genus within Francisellaceae",
            read_ids=["read_010", "read_011"],
        )

    def test_create_valid_cluster(self, novel_species_cluster):
        """Should create NovelCluster from valid data."""
        assert novel_species_cluster.cluster_id == "NSP_001"
        assert novel_species_cluster.taxonomic_call == "Novel Species"
        assert novel_species_cluster.nearest_genome == "GCF_000123456.1"
        assert novel_species_cluster.read_count == 25
        assert novel_species_cluster.confidence == "High"

    def test_cluster_is_frozen(self, novel_species_cluster):
        """NovelCluster should be immutable."""
        with pytest.raises(Exception):
            novel_species_cluster.read_count = 100

    def test_novelty_range_computed_field(self, novel_species_cluster):
        """novelty_range should format min-max correctly."""
        assert novel_species_cluster.novelty_range == "5.2-9.8"

    def test_novelty_band_label_computed_field(self, novel_species_cluster):
        """novelty_band_label should return human-readable label."""
        assert novel_species_cluster.novelty_band_label == "5-10% (Recently diverged)"

    def test_novelty_band_label_moderately_divergent(self, novel_genus_cluster):
        """novelty_band_label for 20% band should indicate genus candidate."""
        assert novel_genus_cluster.novelty_band_label == "20-25% (Novel genus candidate)"

    def test_novelty_band_label_unknown_band(self):
        """novelty_band_label should handle unknown bands."""
        cluster = NovelCluster(
            cluster_id="NSP_999",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_000123456.1",
            nearest_species="Test species",
            nearest_genus="Test",
            nearest_family="Testaceae",
            novelty_band=25,  # Not in predefined labels
            read_count=5,
            mean_novelty_index=27.0,
            novelty_min=25.0,
            novelty_max=29.0,
            mean_placement_uncertainty=2.0,
            suggested_name="Test sp. MDM-999",
            confidence="Low",
            phylogenetic_context="Novel species within Test",
        )
        assert cluster.novelty_band_label == "25-30%"

    def test_to_summary_dict(self, novel_species_cluster):
        """to_summary_dict should return dictionary without read_ids."""
        result = novel_species_cluster.to_summary_dict()

        assert result["cluster_id"] == "NSP_001"
        assert result["taxonomic_call"] == "Novel Species"
        assert result["read_count"] == 25
        assert result["nearest_species"] == "Francisella tularensis"
        assert result["mean_novelty"] == 7.5
        assert result["novelty_range"] == "5.2-9.8"
        assert result["confidence"] == "High"
        assert "read_ids" not in result

    def test_to_summary_dict_rounds_values(self, novel_species_cluster):
        """to_summary_dict should round float values."""
        result = novel_species_cluster.to_summary_dict()

        assert result["mean_uncertainty"] == 1.2
        assert result["mean_discovery_score"] == 82.5

    def test_cluster_without_discovery_score(self):
        """Cluster without discovery_score should work correctly."""
        cluster = NovelCluster(
            cluster_id="NSP_002",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_000123456.1",
            nearest_species="Test species",
            nearest_genus="Test",
            nearest_family="Testaceae",
            novelty_band=10,
            read_count=8,
            mean_novelty_index=12.0,
            novelty_min=10.5,
            novelty_max=14.0,
            mean_placement_uncertainty=4.0,
            mean_discovery_score=None,
            suggested_name="Test sp. MDM-002",
            confidence="Medium",
            phylogenetic_context="Novel species within Test",
        )

        assert cluster.mean_discovery_score is None
        result = cluster.to_summary_dict()
        assert result["mean_discovery_score"] is None

    def test_estimated_ani_computed_field(self, novel_species_cluster):
        """estimated_ani should be 100 - mean_novelty_index."""
        # mean_novelty_index = 7.5, so estimated_ani = 92.5
        assert novel_species_cluster.estimated_ani == 92.5

    def test_closest_known_taxon_novel_species(self, novel_species_cluster):
        """For Novel Species, closest_known_taxon should be nearest_species."""
        assert novel_species_cluster.closest_known_taxon == "Francisella tularensis"

    def test_closest_known_taxon_novel_genus(self, novel_genus_cluster):
        """For Novel Genus, closest_known_taxon should be nearest_genus."""
        assert novel_genus_cluster.closest_known_taxon == "Francisella"

    def test_closest_taxon_with_similarity_novel_species(self, novel_species_cluster):
        """closest_taxon_with_similarity should format correctly for Novel Species."""
        result = novel_species_cluster.closest_taxon_with_similarity
        assert "Francisella tularensis" in result
        assert "92%" in result or "93%" in result  # ~92.5 rounded
        assert "ANI" in result

    def test_closest_taxon_with_similarity_novel_genus(self, novel_genus_cluster):
        """closest_taxon_with_similarity should format correctly for Novel Genus."""
        result = novel_genus_cluster.closest_taxon_with_similarity
        assert "Francisella" in result
        assert "78%" in result  # 100 - 22.3 = 77.7 rounded
        assert "ANI" in result

    def test_to_summary_dict_includes_new_fields(self, novel_species_cluster):
        """to_summary_dict should include closest_known_taxon and estimated_ani."""
        result = novel_species_cluster.to_summary_dict()
        assert "closest_known_taxon" in result
        assert "estimated_ani" in result
        assert result["closest_known_taxon"] == "Francisella tularensis"
        assert result["estimated_ani"] == 92.5


class TestNovelDiversitySummary:
    """Tests for NovelDiversitySummary model."""

    @pytest.fixture
    def sample_summary(self):
        """Sample diversity summary."""
        return NovelDiversitySummary(
            total_novel_reads=200,
            novel_species_reads=150,
            novel_genus_reads=50,
            total_clusters=8,
            novel_species_clusters=5,
            novel_genus_clusters=3,
            high_confidence_clusters=2,
            medium_confidence_clusters=4,
            low_confidence_clusters=2,
            genera_with_novel_species=3,
            families_with_novel_genera=2,
            mean_cluster_size=25.0,
            largest_cluster_size=50,
        )

    def test_create_valid_summary(self, sample_summary):
        """Should create NovelDiversitySummary from valid data."""
        assert sample_summary.total_novel_reads == 200
        assert sample_summary.total_clusters == 8
        assert sample_summary.high_confidence_clusters == 2

    def test_default_values(self):
        """Summary should have sensible defaults."""
        summary = NovelDiversitySummary()

        assert summary.total_novel_reads == 0
        assert summary.total_clusters == 0
        assert summary.mean_cluster_size == 0.0

    def test_novel_species_pct_computed_field(self, sample_summary):
        """novel_species_pct should calculate correctly."""
        # 150 / 200 * 100 = 75%
        assert sample_summary.novel_species_pct == 75.0

    def test_novel_genus_pct_computed_field(self, sample_summary):
        """novel_genus_pct should calculate correctly."""
        # 50 / 200 * 100 = 25%
        assert sample_summary.novel_genus_pct == 25.0

    def test_percentages_with_zero_reads(self):
        """Percentages should be 0 when total_novel_reads is 0."""
        summary = NovelDiversitySummary()

        assert summary.novel_species_pct == 0.0
        assert summary.novel_genus_pct == 0.0


# =============================================================================
# NovelDiversityAnalyzer Tests
# =============================================================================


class TestNovelDiversityAnalyzer:
    """Tests for NovelDiversityAnalyzer clustering logic."""

    @pytest.fixture
    def classification_df_with_novel(self):
        """Classification DataFrame with novel reads."""
        return pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(20)],
            "best_match_genome": (
                ["GCF_000123456.1"] * 10 +
                ["GCF_000789012.1"] * 6 +
                ["GCA_000111222.1"] * 4
            ),
            "novelty_index": (
                [6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 12.0, 13.0, 14.0] +  # Novel Species
                [21.0, 22.0, 23.0, 24.0, 8.5, 9.0] +  # Novel Genus + Novel Species
                [1.0, 2.0, 3.0, 15.0]  # Known Species + Conserved
            ),
            "placement_uncertainty": (
                [1.0] * 10 + [2.0] * 6 + [0.5, 0.6, 0.7, 12.0]
            ),
            "taxonomic_call": (
                ["Novel Species"] * 10 +
                ["Novel Genus"] * 4 + ["Novel Species"] * 2 +
                ["Known Species"] * 3 + ["Conserved Region"]
            ),
        })

    @pytest.fixture
    def classification_df_no_novel(self):
        """Classification DataFrame without novel reads."""
        return pl.DataFrame({
            "read_id": ["read_001", "read_002", "read_003"],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [1.0, 2.0, 3.0],
            "placement_uncertainty": [0.5, 0.6, 0.7],
            "taxonomic_call": ["Known Species"] * 3,
        })

    @pytest.fixture
    def classification_df_with_discovery_score(self):
        """Classification DataFrame with discovery_score column."""
        return pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(10)],
            "best_match_genome": ["GCF_000123456.1"] * 10,
            "novelty_index": [7.0 + i * 0.5 for i in range(10)],
            "placement_uncertainty": [3.0] * 10,
            "taxonomic_call": ["Novel Species"] * 10,
            "discovery_score": [80.0, 82.0, 78.0, 85.0, 90.0, 75.0, 88.0, 92.0, 70.0, 77.0],
        })

    def test_init_validates_required_columns(self):
        """Analyzer should validate required columns."""
        df_missing_cols = pl.DataFrame({
            "read_id": ["read_001"],
            "best_match_genome": ["GCF_000123456.1"],
        })

        with pytest.raises(ValueError, match="missing required columns"):
            NovelDiversityAnalyzer(df_missing_cols)

    def test_cluster_novel_reads_basic(self, classification_df_with_novel):
        """Should cluster novel reads into groups."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=3,
        )

        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) >= 1
        # All clusters should have at least min_cluster_size reads
        for cluster in clusters:
            assert cluster.read_count >= 3

    def test_cluster_novel_reads_no_novel(self, classification_df_no_novel):
        """Should return empty list when no novel reads exist."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_no_novel,
            min_cluster_size=1,
        )

        clusters = analyzer.cluster_novel_reads()

        assert clusters == []

    def test_cluster_caching(self, classification_df_with_novel):
        """Clustering should be cached after first call."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=3,
        )

        clusters1 = analyzer.cluster_novel_reads()
        clusters2 = analyzer.cluster_novel_reads()

        assert clusters1 is clusters2

    def test_cluster_ids_are_unique(self, classification_df_with_novel):
        """All cluster IDs should be unique."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=1,
        )

        clusters = analyzer.cluster_novel_reads()
        cluster_ids = [c.cluster_id for c in clusters]

        assert len(cluster_ids) == len(set(cluster_ids))

    def test_cluster_id_format(self, classification_df_with_novel):
        """Cluster IDs should follow NSP_XXX or NGN_XXX format."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=1,
        )

        clusters = analyzer.cluster_novel_reads()

        for cluster in clusters:
            if cluster.taxonomic_call == "Novel Species":
                assert cluster.cluster_id.startswith("NSP_")
            else:
                assert cluster.cluster_id.startswith("NGN_")
            # Should have 3-digit suffix
            suffix = cluster.cluster_id.split("_")[1]
            assert len(suffix) == 3
            assert suffix.isdigit()

    def test_novelty_band_calculation(self, classification_df_with_novel):
        """Novelty bands should group reads correctly."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            novelty_band_size=5.0,
            min_cluster_size=1,
        )

        clusters = analyzer.cluster_novel_reads()

        for cluster in clusters:
            # Band should be multiple of 5
            assert cluster.novelty_band % 5 == 0
            # Mean novelty should fall within band range
            assert cluster.novelty_band <= cluster.mean_novelty_index
            assert cluster.mean_novelty_index < cluster.novelty_band + 5

    def test_discovery_score_aggregation(self, classification_df_with_discovery_score):
        """Should aggregate discovery_score when present."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_discovery_score,
            min_cluster_size=1,
        )

        clusters = analyzer.cluster_novel_reads()

        for cluster in clusters:
            assert cluster.mean_discovery_score is not None
            assert 0 <= cluster.mean_discovery_score <= 100

    def test_get_summary(self, classification_df_with_novel):
        """get_summary should return aggregate statistics."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=3,
        )

        summary = analyzer.get_summary()

        assert isinstance(summary, NovelDiversitySummary)
        assert summary.total_novel_reads > 0
        # Novel species reads + novel genus reads = total
        assert (
            summary.novel_species_reads + summary.novel_genus_reads
            == summary.total_novel_reads
        )

    def test_get_summary_caching(self, classification_df_with_novel):
        """Summary should be cached after first call."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=3,
        )

        summary1 = analyzer.get_summary()
        summary2 = analyzer.get_summary()

        assert summary1 is summary2

    def test_get_summary_empty(self, classification_df_no_novel):
        """Summary should have zeros when no novel reads."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_no_novel,
            min_cluster_size=1,
        )

        summary = analyzer.get_summary()

        assert summary.total_novel_reads == 0
        assert summary.total_clusters == 0

    def test_to_dataframe(self, classification_df_with_novel):
        """to_dataframe should export clusters as DataFrame."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=3,
        )

        df = analyzer.to_dataframe()

        assert isinstance(df, pl.DataFrame)
        if len(df) > 0:
            assert "cluster_id" in df.columns
            assert "taxonomic_call" in df.columns
            assert "read_count" in df.columns

    def test_to_dataframe_empty(self, classification_df_no_novel):
        """to_dataframe should return empty DataFrame when no clusters."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_no_novel,
            min_cluster_size=1,
        )

        df = analyzer.to_dataframe()

        assert isinstance(df, pl.DataFrame)
        assert len(df) == 0

    def test_to_dict(self, classification_df_with_novel):
        """to_dict should export as JSON-serializable dictionary."""
        analyzer = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=3,
        )

        result = analyzer.to_dict()

        assert "summary" in result
        assert "clusters" in result
        assert isinstance(result["summary"], dict)
        assert isinstance(result["clusters"], list)

    def test_min_cluster_size_filtering(self, classification_df_with_novel):
        """min_cluster_size should filter out small clusters."""
        analyzer_low = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=1,
        )
        analyzer_high = NovelDiversityAnalyzer(
            classifications=classification_df_with_novel,
            min_cluster_size=10,
        )

        clusters_low = analyzer_low.cluster_novel_reads()
        clusters_high = analyzer_high.cluster_novel_reads()

        # Higher min_cluster_size should result in fewer or equal clusters
        assert len(clusters_high) <= len(clusters_low)


class TestConfidenceRating:
    """Tests for confidence rating assignment."""

    @pytest.fixture
    def analyzer_with_discovery(self):
        """Analyzer with discovery scores for confidence testing."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(30)],
            "best_match_genome": ["GCF_000123456.1"] * 30,
            "novelty_index": [7.0] * 30,
            "placement_uncertainty": [2.0] * 10 + [8.0] * 10 + [4.0] * 10,
            "taxonomic_call": ["Novel Species"] * 30,
            "discovery_score": [80.0] * 10 + [60.0] * 10 + [40.0] * 10,
        })
        return NovelDiversityAnalyzer(df, min_cluster_size=1)

    def test_high_confidence_criteria(self):
        """High confidence requires high read count, low uncertainty, high discovery."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(15)],
            "best_match_genome": ["GCF_000123456.1"] * 15,
            "novelty_index": [7.0] * 15,
            "placement_uncertainty": [2.0] * 15,  # < 5%
            "taxonomic_call": ["Novel Species"] * 15,
            "discovery_score": [80.0] * 15,  # >= 75
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=10)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "High"

    def test_medium_confidence_criteria(self):
        """Medium confidence with moderate criteria."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(7)],
            "best_match_genome": ["GCF_000123456.1"] * 7,
            "novelty_index": [7.0] * 7,
            "placement_uncertainty": [7.0] * 7,  # < 10%
            "taxonomic_call": ["Novel Species"] * 7,
            "discovery_score": [55.0] * 7,  # >= 50
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=5)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "Medium"

    def test_low_confidence_criteria(self):
        """Low confidence for clusters not meeting other criteria."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(3)],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [7.0] * 3,
            "placement_uncertainty": [15.0] * 3,  # High uncertainty
            "taxonomic_call": ["Novel Species"] * 3,
            "discovery_score": [30.0] * 3,  # Low discovery
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=3)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "Low"

    def test_confidence_without_discovery_score(self):
        """Confidence should work without discovery_score."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(15)],
            "best_match_genome": ["GCF_000123456.1"] * 15,
            "novelty_index": [7.0] * 15,
            "placement_uncertainty": [2.0] * 15,  # < 5%
            "taxonomic_call": ["Novel Species"] * 15,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=10)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "High"


class TestPhylogeneticContext:
    """Tests for phylogenetic context inference."""

    @pytest.fixture
    def mock_metadata(self):
        """Mock genome metadata."""
        from unittest.mock import Mock

        metadata = Mock()
        metadata.get_species.return_value = "Francisella tularensis"
        metadata.get_genus.return_value = "Francisella"
        metadata.get_family.return_value = "Francisellaceae"
        return metadata

    def test_novel_species_context_with_genus(self, mock_metadata):
        """Novel Species should show 'within {genus}'."""
        df = pl.DataFrame({
            "read_id": ["read_001", "read_002", "read_003"],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [7.0] * 3,
            "placement_uncertainty": [2.0] * 3,
            "taxonomic_call": ["Novel Species"] * 3,
        })
        analyzer = NovelDiversityAnalyzer(
            df,
            metadata=mock_metadata,
            min_cluster_size=1,
        )
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert "Novel species within Francisella" in clusters[0].phylogenetic_context

    def test_novel_genus_context_with_family(self, mock_metadata):
        """Novel Genus should show 'within {family}'."""
        df = pl.DataFrame({
            "read_id": ["read_001", "read_002", "read_003"],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [22.0] * 3,
            "placement_uncertainty": [2.0] * 3,
            "taxonomic_call": ["Novel Genus"] * 3,
        })
        analyzer = NovelDiversityAnalyzer(
            df,
            metadata=mock_metadata,
            min_cluster_size=1,
        )
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert "Novel genus within Francisellaceae" in clusters[0].phylogenetic_context

    def test_context_without_metadata(self):
        """Context should work without metadata (using Unknown)."""
        df = pl.DataFrame({
            "read_id": ["read_001", "read_002", "read_003"],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [7.0] * 3,
            "placement_uncertainty": [2.0] * 3,
            "taxonomic_call": ["Novel Species"] * 3,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=1)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert "genus unknown" in clusters[0].phylogenetic_context.lower()


class TestInferredUncertainty:
    """Tests for novelty-based uncertainty inference for single hits."""

    def test_single_hit_uncertainty_inferred(self):
        """Single-hit reads should have uncertainty inferred from novelty."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(5)],
            "best_match_genome": ["GCF_000123456.1"] * 5,
            "novelty_index": [10.0] * 5,  # 10% novelty
            "placement_uncertainty": [0.0] * 5,  # Single hits have 0 uncertainty
            "num_ambiguous_hits": [1] * 5,  # Single hits
            "taxonomic_call": ["Novel Species"] * 5,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=3)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        # Inferred uncertainty should be novelty * 0.4 = 10 * 0.4 = 4.0
        assert clusters[0].mean_placement_uncertainty == pytest.approx(4.0, rel=0.1)

    def test_multi_hit_uncertainty_preserved(self):
        """Multi-hit reads should keep their actual uncertainty."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(5)],
            "best_match_genome": ["GCF_000123456.1"] * 5,
            "novelty_index": [10.0] * 5,
            "placement_uncertainty": [3.5] * 5,  # Actual uncertainty from multiple hits
            "num_ambiguous_hits": [5] * 5,  # Multiple hits
            "taxonomic_call": ["Novel Species"] * 5,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=3)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        # Should use actual uncertainty, not inferred
        assert clusters[0].mean_placement_uncertainty == pytest.approx(3.5, rel=0.1)

    def test_inferred_uncertainty_capped(self):
        """Inferred uncertainty should be capped at 15%."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(5)],
            "best_match_genome": ["GCF_000123456.1"] * 5,
            "novelty_index": [50.0] * 5,  # High novelty would give 20% inferred
            "placement_uncertainty": [0.0] * 5,
            "num_ambiguous_hits": [1] * 5,
            "taxonomic_call": ["Novel Genus"] * 5,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=3)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        # Should be capped at 15%, not 50 * 0.4 = 20
        assert clusters[0].mean_placement_uncertainty == pytest.approx(15.0, rel=0.1)

    def test_fallback_without_num_ambiguous_hits(self):
        """Should infer for zero-uncertainty reads when num_ambiguous_hits missing."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(5)],
            "best_match_genome": ["GCF_000123456.1"] * 5,
            "novelty_index": [8.0] * 5,
            "placement_uncertainty": [0.0] * 5,  # Zero uncertainty triggers inference
            "taxonomic_call": ["Novel Species"] * 5,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=3)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        # Inferred: 8.0 * 0.4 = 3.2
        assert clusters[0].mean_placement_uncertainty == pytest.approx(3.2, rel=0.1)


class TestSuggestedNames:
    """Tests for GTDB-style provisional name generation."""

    @pytest.fixture
    def mock_metadata(self):
        """Mock genome metadata."""
        from unittest.mock import Mock

        metadata = Mock()
        metadata.get_species.return_value = "Francisella tularensis"
        metadata.get_genus.return_value = "Francisella"
        metadata.get_family.return_value = "Francisellaceae"
        return metadata

    def test_novel_species_name_format(self, mock_metadata):
        """Novel Species name should follow '{Genus} sp. MDM-XXX' format."""
        df = pl.DataFrame({
            "read_id": ["read_001", "read_002", "read_003"],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [7.0] * 3,
            "placement_uncertainty": [2.0] * 3,
            "taxonomic_call": ["Novel Species"] * 3,
        })
        analyzer = NovelDiversityAnalyzer(
            df,
            metadata=mock_metadata,
            min_cluster_size=1,
        )
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].suggested_name.startswith("Francisella sp. MDM-")

    def test_novel_genus_name_format_with_family(self, mock_metadata):
        """Novel Genus name should follow '{Family} gen. nov. MDM-XXX' format."""
        df = pl.DataFrame({
            "read_id": ["read_001", "read_002", "read_003"],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [22.0] * 3,
            "placement_uncertainty": [2.0] * 3,
            "taxonomic_call": ["Novel Genus"] * 3,
        })
        analyzer = NovelDiversityAnalyzer(
            df,
            metadata=mock_metadata,
            min_cluster_size=1,
        )
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert "Francisellaceae gen. nov. MDM-" in clusters[0].suggested_name

    def test_name_without_metadata(self):
        """Name should use Unknown genus/family without metadata."""
        df = pl.DataFrame({
            "read_id": ["read_001", "read_002", "read_003"],
            "best_match_genome": ["GCF_000123456.1"] * 3,
            "novelty_index": [7.0] * 3,
            "placement_uncertainty": [2.0] * 3,
            "taxonomic_call": ["Novel Species"] * 3,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=1)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].suggested_name.startswith("Unknown sp. MDM-")


# =============================================================================
# NovelCluster New Fields Tests
# =============================================================================


class TestNovelClusterNewFields:
    """Tests for mean_bayesian_confidence and contributing_genomes fields."""

    def test_new_fields_default_to_none_and_empty(self):
        """New fields should default gracefully."""
        cluster = NovelCluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_000123456.1",
            nearest_species="Test species",
            nearest_genus="Test",
            nearest_family="Testaceae",
            novelty_band=5,
            read_count=10,
            mean_novelty_index=7.0,
            novelty_min=5.0,
            novelty_max=9.0,
            mean_placement_uncertainty=2.0,
            suggested_name="Test sp. MDM-001",
            confidence="High",
            phylogenetic_context="Novel species within Test",
        )
        assert cluster.mean_bayesian_confidence is None
        assert cluster.contributing_genomes == []

    def test_new_fields_can_be_set(self):
        """New fields should accept values."""
        cluster = NovelCluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_000123456.1",
            nearest_species="Test species",
            nearest_genus="Test",
            nearest_family="Testaceae",
            novelty_band=5,
            read_count=10,
            mean_novelty_index=7.0,
            novelty_min=5.0,
            novelty_max=9.0,
            mean_placement_uncertainty=2.0,
            suggested_name="Test sp. MDM-001",
            confidence="High",
            phylogenetic_context="Novel species within Test",
            mean_bayesian_confidence=72.5,
            contributing_genomes=["GCF_A", "GCF_B"],
        )
        assert cluster.mean_bayesian_confidence == 72.5
        assert cluster.contributing_genomes == ["GCF_A", "GCF_B"]

    def test_to_summary_dict_includes_bayesian_confidence(self):
        """to_summary_dict should include mean_bayesian_confidence when set."""
        cluster = NovelCluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_000123456.1",
            nearest_species="Test species",
            nearest_genus="Test",
            nearest_family="Testaceae",
            novelty_band=5,
            read_count=10,
            mean_novelty_index=7.0,
            novelty_min=5.0,
            novelty_max=9.0,
            mean_placement_uncertainty=2.0,
            suggested_name="Test sp. MDM-001",
            confidence="High",
            phylogenetic_context="Novel species within Test",
            mean_bayesian_confidence=72.5,
            contributing_genomes=["GCF_A", "GCF_B"],
        )
        result = cluster.to_summary_dict()
        assert result["mean_bayesian_confidence"] == 72.5
        assert result["contributing_genomes_count"] == 2

    def test_to_summary_dict_omits_unset_new_fields(self):
        """to_summary_dict should omit new fields when not set."""
        cluster = NovelCluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_000123456.1",
            nearest_species="Test species",
            nearest_genus="Test",
            nearest_family="Testaceae",
            novelty_band=5,
            read_count=10,
            mean_novelty_index=7.0,
            novelty_min=5.0,
            novelty_max=9.0,
            mean_placement_uncertainty=2.0,
            suggested_name="Test sp. MDM-001",
            confidence="High",
            phylogenetic_context="Novel species within Test",
        )
        result = cluster.to_summary_dict()
        assert "mean_bayesian_confidence" not in result
        assert "contributing_genomes_count" not in result


# =============================================================================
# Genome Neighborhood Merging Tests
# =============================================================================


class TestGenomeNeighborhoods:
    """Tests for ANI-aware genome neighborhood merging."""

    @pytest.fixture
    def ani_matrix(self):
        """ANI matrix with two neighborhoods."""
        from metadarkmatter.core.classification.ani_matrix import ANIMatrix

        ani_dict = {
            "GCF_A": {"GCF_A": 100.0, "GCF_B": 87.0, "GCF_C": 72.0},
            "GCF_B": {"GCF_A": 87.0, "GCF_B": 100.0, "GCF_C": 73.0},
            "GCF_C": {"GCF_A": 72.0, "GCF_B": 73.0, "GCF_C": 100.0},
        }
        return ANIMatrix(ani_dict)

    def test_neighborhood_merging_reduces_clusters(self, ani_matrix):
        """Genomes within ANI threshold should merge into one neighborhood."""
        # GCF_A and GCF_B are at 87% ANI (above 80% threshold)
        # So reads hitting either should cluster together
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(20)],
            "best_match_genome": ["GCF_A"] * 10 + ["GCF_B"] * 10,
            "novelty_index": [7.0] * 20,
            "placement_uncertainty": [2.0] * 20,
            "taxonomic_call": ["Novel Species"] * 20,
        })

        # Without ANI: 2 clusters (one per genome)
        analyzer_no_ani = NovelDiversityAnalyzer(
            df, min_cluster_size=3,
        )
        clusters_no_ani = analyzer_no_ani.cluster_novel_reads()

        # With ANI: 1 cluster (merged neighborhood)
        analyzer_with_ani = NovelDiversityAnalyzer(
            df, ani_matrix=ani_matrix, min_cluster_size=3,
        )
        clusters_with_ani = analyzer_with_ani.cluster_novel_reads()

        assert len(clusters_no_ani) == 2
        assert len(clusters_with_ani) == 1
        assert clusters_with_ani[0].read_count == 20

    def test_separate_neighborhoods_stay_separate(self, ani_matrix):
        """Genomes below ANI threshold should remain in separate clusters."""
        # GCF_A and GCF_C are at 72% ANI (below 80% threshold)
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(10)],
            "best_match_genome": ["GCF_A"] * 5 + ["GCF_C"] * 5,
            "novelty_index": [7.0] * 10,
            "placement_uncertainty": [2.0] * 10,
            "taxonomic_call": ["Novel Species"] * 10,
        })

        analyzer = NovelDiversityAnalyzer(
            df, ani_matrix=ani_matrix, min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 2

    def test_contributing_genomes_populated(self, ani_matrix):
        """contributing_genomes should list all original genomes in merged cluster."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(20)],
            "best_match_genome": ["GCF_A"] * 10 + ["GCF_B"] * 10,
            "novelty_index": [7.0] * 20,
            "placement_uncertainty": [2.0] * 20,
            "taxonomic_call": ["Novel Species"] * 20,
        })

        analyzer = NovelDiversityAnalyzer(
            df, ani_matrix=ani_matrix, min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert set(clusters[0].contributing_genomes) == {"GCF_A", "GCF_B"}

    def test_no_ani_matrix_preserves_old_behavior(self):
        """Without ANI matrix, clustering should work as before."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(10)],
            "best_match_genome": ["GCF_A"] * 5 + ["GCF_B"] * 5,
            "novelty_index": [7.0] * 10,
            "placement_uncertainty": [2.0] * 10,
            "taxonomic_call": ["Novel Species"] * 10,
        })

        analyzer = NovelDiversityAnalyzer(
            df, min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        # Two genomes should produce two clusters without ANI merging
        assert len(clusters) == 2


# =============================================================================
# Adjacent Band Merging Tests
# =============================================================================


class TestAdjacentBandMerging:
    """Tests for adjacent novelty band merging."""

    def test_adjacent_bands_merge(self):
        """Reads at 9.9% and 10.1% novelty should merge into one cluster."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(20)],
            "best_match_genome": ["GCF_A"] * 20,
            "novelty_index": [9.0, 9.5, 9.8, 9.9, 9.9] * 2 + [10.1, 10.2, 10.5, 11.0, 11.5] * 2,
            "placement_uncertainty": [2.0] * 20,
            "taxonomic_call": ["Novel Species"] * 20,
        })

        analyzer = NovelDiversityAnalyzer(
            df, min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        # Should merge into 1 cluster (bands 5 and 10 are adjacent, ranges overlap)
        assert len(clusters) == 1
        assert clusters[0].read_count == 20

    def test_non_adjacent_bands_stay_separate(self):
        """Reads in bands 5 and 15 should not merge (not adjacent)."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(20)],
            "best_match_genome": ["GCF_A"] * 20,
            "novelty_index": [6.0] * 10 + [16.0] * 10,
            "placement_uncertainty": [2.0] * 20,
            "taxonomic_call": ["Novel Species"] * 20,
        })

        analyzer = NovelDiversityAnalyzer(
            df, min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        # Should remain 2 separate clusters (bands 5 and 15 are not adjacent)
        assert len(clusters) == 2

    def test_different_tax_calls_dont_merge(self):
        """Adjacent bands with different taxonomic calls should not merge."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(10)],
            "best_match_genome": ["GCF_A"] * 10,
            "novelty_index": [9.0] * 5 + [10.5] * 5,
            "placement_uncertainty": [2.0] * 10,
            "taxonomic_call": ["Novel Species"] * 5 + ["Novel Genus"] * 5,
        })

        analyzer = NovelDiversityAnalyzer(
            df, min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        # Different taxonomic calls: should not merge
        assert len(clusters) == 2


# =============================================================================
# Bayesian Confidence Tests
# =============================================================================


class TestBayesianConfidence:
    """Tests for Bayesian-informed confidence assignment."""

    def test_high_bayesian_confidence(self):
        """High Bayesian confidence (>=70) with >=10 reads should be High."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(15)],
            "best_match_genome": ["GCF_A"] * 15,
            "novelty_index": [7.0] * 15,
            "placement_uncertainty": [2.0] * 15,
            "taxonomic_call": ["Novel Species"] * 15,
            "confidence_score": [75.0] * 15,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=10)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "High"
        assert clusters[0].mean_bayesian_confidence == pytest.approx(75.0, rel=0.1)

    def test_medium_bayesian_confidence(self):
        """Bayesian confidence 50-70 with >=5 reads should be Medium."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(7)],
            "best_match_genome": ["GCF_A"] * 7,
            "novelty_index": [7.0] * 7,
            "placement_uncertainty": [2.0] * 7,
            "taxonomic_call": ["Novel Species"] * 7,
            "confidence_score": [55.0] * 7,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=5)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "Medium"

    def test_low_bayesian_confidence(self):
        """Low Bayesian confidence (<50) should be Low regardless of reads."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(20)],
            "best_match_genome": ["GCF_A"] * 20,
            "novelty_index": [7.0] * 20,
            "placement_uncertainty": [2.0] * 20,
            "taxonomic_call": ["Novel Species"] * 20,
            "confidence_score": [30.0] * 20,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=10)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "Low"

    def test_bayesian_overrides_legacy(self):
        """Bayesian confidence should override legacy discovery-based criteria."""
        # Legacy criteria would give High (>=10 reads, uncertainty <5%, discovery >=75)
        # But Bayesian confidence is 40 -> should be Low
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(15)],
            "best_match_genome": ["GCF_A"] * 15,
            "novelty_index": [7.0] * 15,
            "placement_uncertainty": [2.0] * 15,
            "taxonomic_call": ["Novel Species"] * 15,
            "discovery_score": [80.0] * 15,
            "confidence_score": [40.0] * 15,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=10)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        assert clusters[0].confidence == "Low"

    def test_fallback_without_confidence_score(self):
        """Without confidence_score column, should use legacy criteria."""
        df = pl.DataFrame({
            "read_id": [f"read_{i:03d}" for i in range(15)],
            "best_match_genome": ["GCF_A"] * 15,
            "novelty_index": [7.0] * 15,
            "placement_uncertainty": [2.0] * 15,
            "taxonomic_call": ["Novel Species"] * 15,
        })
        analyzer = NovelDiversityAnalyzer(df, min_cluster_size=10)
        clusters = analyzer.cluster_novel_reads()

        assert len(clusters) == 1
        # Legacy: >=10 reads, uncertainty < 5% -> High
        assert clusters[0].confidence == "High"
        assert clusters[0].mean_bayesian_confidence is None


# =============================================================================
# ANIMatrix.from_dataframe Tests
# =============================================================================


class TestANIMatrixFromDataframe:
    """Tests for ANIMatrix.from_dataframe() classmethod."""

    def test_from_dataframe_with_genome_column(self):
        """Should construct from DataFrame with 'genome' column."""
        import polars as pl
        from metadarkmatter.core.classification.ani_matrix import ANIMatrix

        df = pl.DataFrame({
            "genome": ["GCF_A", "GCF_B", "GCF_C"],
            "GCF_A": [100.0, 96.0, 82.0],
            "GCF_B": [96.0, 100.0, 83.0],
            "GCF_C": [82.0, 83.0, 100.0],
        })

        matrix = ANIMatrix.from_dataframe(df)

        assert len(matrix) == 3
        assert matrix.get_ani("GCF_A", "GCF_B") == pytest.approx(96.0)
        assert matrix.get_ani("GCF_A", "GCF_C") == pytest.approx(82.0)
        assert matrix.get_ani("GCF_A", "GCF_A") == 100.0

    def test_from_dataframe_without_genome_column(self):
        """Should construct from square DataFrame without 'genome' column."""
        import polars as pl
        from metadarkmatter.core.classification.ani_matrix import ANIMatrix

        df = pl.DataFrame({
            "GCF_A": [100.0, 96.0],
            "GCF_B": [96.0, 100.0],
        })

        matrix = ANIMatrix.from_dataframe(df)

        assert len(matrix) == 2
        assert matrix.get_ani("GCF_A", "GCF_B") == pytest.approx(96.0)

    def test_from_dataframe_empty(self):
        """Should handle empty DataFrame."""
        import polars as pl
        from metadarkmatter.core.classification.ani_matrix import ANIMatrix

        df = pl.DataFrame()

        matrix = ANIMatrix.from_dataframe(df)

        assert len(matrix) == 0

    def test_from_dataframe_roundtrip(self):
        """Values should survive dict -> ANIMatrix -> from_dataframe cycle."""
        import polars as pl
        from metadarkmatter.core.classification.ani_matrix import ANIMatrix

        ani_dict = {
            "GCF_A": {"GCF_A": 100.0, "GCF_B": 95.5, "GCF_C": 80.0},
            "GCF_B": {"GCF_A": 95.5, "GCF_B": 100.0, "GCF_C": 81.0},
            "GCF_C": {"GCF_A": 80.0, "GCF_B": 81.0, "GCF_C": 100.0},
        }
        original = ANIMatrix(ani_dict)

        # Build a DataFrame matching the report generator's format
        genomes = sorted(ani_dict.keys())
        data = {"genome": genomes}
        for col_genome in genomes:
            data[col_genome] = [original.get_ani(row_g, col_genome) for row_g in genomes]
        df = pl.DataFrame(data)

        reconstructed = ANIMatrix.from_dataframe(df)

        for g1 in genomes:
            for g2 in genomes:
                assert reconstructed.get_ani(g1, g2) == pytest.approx(
                    original.get_ani(g1, g2), abs=0.1
                )
