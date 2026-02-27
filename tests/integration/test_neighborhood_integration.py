"""Integration tests for the phylogenetic neighborhood pipeline.

Exercises the full pipeline from ANI matrix construction through
classification, clustering, neighborhood analysis, and report generation
to verify that all stages connect correctly and produce expected output.
"""

from __future__ import annotations

import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.novel_diversity import (
    NovelDiversityAnalyzer,
    PhylogeneticNeighborhoodAnalyzer,
)

# Skip the module when plotly is absent (needed for report generation tests).
plotly = pytest.importorskip("plotly")


# =============================================================================
# Helpers
# =============================================================================


def _make_multi_genus_ani_matrix() -> tuple[ANIMatrix, list[str]]:
    """Create a 9-genome ANI matrix with three genera.

    Genus layout:
        GenusA: GCF_001, GCF_002, GCF_003  (within-genus ~92% ANI)
        GenusB: GCF_004, GCF_005, GCF_006  (within-genus ~92% ANI)
        GenusC: GCF_007, GCF_008, GCF_009  (within-genus ~92% ANI)
    Inter-genus ANI is fixed at 75%.

    Returns:
        Tuple of (ANIMatrix, genome_list).
    """
    genomes = [f"GCF_{i:03d}" for i in range(1, 10)]
    n = len(genomes)

    ani_dict: dict[str, dict[str, float]] = {}
    for i, g1 in enumerate(genomes):
        inner: dict[str, float] = {}
        for j, g2 in enumerate(genomes):
            if i == j:
                inner[g2] = 100.0
            elif i // 3 == j // 3:
                inner[g2] = 92.0
            else:
                inner[g2] = 75.0
        ani_dict[g1] = inner

    return ANIMatrix(ani_dict), genomes


def _make_genus_map(genomes: list[str]) -> dict[str, str]:
    """Map genomes to three genera based on index."""
    genera = ["GenusA", "GenusB", "GenusC"]
    return {g: genera[i // 3] for i, g in enumerate(genomes)}


def _make_novel_classifications(genomes: list[str]) -> pl.DataFrame:
    """Create a classification DataFrame with novel and known reads.

    Generates:
        - 20 Novel Genus reads hitting GCF_001 (GenusA)
        - 15 Novel Species reads hitting GCF_004 (GenusB)
        - 30 Known Species reads hitting GCF_007 (GenusC)

    Args:
        genomes: Ordered genome accession list from _make_multi_genus_ani_matrix.

    Returns:
        Polars DataFrame matching the schema expected by NovelDiversityAnalyzer.
    """
    rows: list[dict] = []

    # Novel Genus reads -- high novelty, hitting GenusA representative
    for i in range(20):
        rows.append({
            "read_id": f"novel_genus_{i}",
            "best_match_genome": genomes[0],  # GCF_001 (GenusA)
            "top_hit_identity": 78.0,
            "novelty_index": 22.0,
            "placement_uncertainty": 1.0,
            "taxonomic_call": "Novel Genus",
            "num_ambiguous_hits": 1,
        })

    # Novel Species reads -- moderate novelty, hitting GenusB representative
    for i in range(15):
        rows.append({
            "read_id": f"novel_species_{i}",
            "best_match_genome": genomes[3],  # GCF_004 (GenusB)
            "top_hit_identity": 91.0,
            "novelty_index": 9.0,
            "placement_uncertainty": 0.8,
            "taxonomic_call": "Novel Species",
            "num_ambiguous_hits": 1,
        })

    # Known Species reads -- low novelty
    for i in range(30):
        rows.append({
            "read_id": f"known_{i}",
            "best_match_genome": genomes[6],  # GCF_007 (GenusC)
            "top_hit_identity": 98.0,
            "novelty_index": 2.0,
            "placement_uncertainty": 0.5,
            "taxonomic_call": "Known Species",
            "num_ambiguous_hits": 2,
        })

    return pl.DataFrame(rows)


def _make_metadata_polars(genomes: list[str]) -> pl.DataFrame:
    """Build a metadata Polars DataFrame suitable for the ReportGenerator.

    Provides genome, genus, and species columns for report-level operations.
    """
    genera = ["GenusA", "GenusB", "GenusC"]
    rows = []
    for i, g in enumerate(genomes):
        genus = genera[i // 3]
        rows.append({
            "accession": g,
            "species": f"{genus} sp{i % 3 + 1}",
            "genus": genus,
            "family": "TestFamily",
        })
    return pl.DataFrame(rows)


# =============================================================================
# Integration tests
# =============================================================================


class TestNeighborhoodIntegration:
    """Full pipeline: ANI matrix -> clustering -> neighborhood."""

    def test_full_pipeline_clustering_to_neighborhood(self):
        """Clustering followed by neighborhood analysis should produce
        enriched clusters with populated neighborhood fields."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        # Step 1: Cluster novel reads
        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()
        assert len(clusters) > 0

        # Step 2: Neighborhood enrichment
        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
            genus_boundary=80.0,
        )
        enriched = nbr_analyzer.analyze(clusters)

        # Verify every cluster gained a neighborhood
        for cluster in enriched:
            assert cluster.neighborhood is not None
            assert len(cluster.neighborhood.nearest_genera) >= 1
            assert 0 <= cluster.neighborhood.placement_support <= 100
            assert cluster.neighborhood.phylogenetic_context  # non-empty

    def test_novel_genus_nearest_genus_is_correct(self):
        """Novel Genus cluster should report GenusA as its nearest genus,
        since its reads hit GCF_001 which belongs to GenusA."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        enriched = nbr_analyzer.analyze(clusters)

        genus_clusters = [
            c for c in enriched if c.taxonomic_call == "Novel Genus"
        ]
        assert len(genus_clusters) >= 1, "Expected at least one Novel Genus cluster"

        nbr = genus_clusters[0].neighborhood
        assert nbr is not None
        assert nbr.nearest_genera[0].genus == "GenusA"

    def test_novel_species_nearest_genus_is_correct(self):
        """Novel Species cluster should report GenusB as nearest genus,
        since its reads hit GCF_004 which belongs to GenusB."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        enriched = nbr_analyzer.analyze(clusters)

        species_clusters = [
            c for c in enriched if c.taxonomic_call == "Novel Species"
        ]
        assert len(species_clusters) >= 1, "Expected at least one Novel Species cluster"

        nbr = species_clusters[0].neighborhood
        assert nbr is not None
        assert nbr.nearest_genera[0].genus == "GenusB"

    def test_context_text_reflects_taxonomy(self):
        """Context text should reference ANI and placement support."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        enriched = nbr_analyzer.analyze(clusters)

        for cluster in enriched:
            nbr = cluster.neighborhood
            if nbr is not None:
                assert "ANI" in nbr.phylogenetic_context
                assert "Support" in nbr.phylogenetic_context

    def test_isolation_and_density_metrics(self):
        """Isolation score and neighborhood density should be consistent
        with the three-genera reference layout."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
            genus_boundary=80.0,
        )
        enriched = nbr_analyzer.analyze(clusters)

        for cluster in enriched:
            nbr = cluster.neighborhood
            assert nbr is not None
            # With 3 genera all at known ANI values, isolation score should be >= 0
            assert nbr.isolation_score >= 0.0
            # Density: at least 1 genus within 5% of genus boundary
            assert nbr.neighborhood_density >= 1


class TestNeighborhoodReportIntegration:
    """Pipeline through to report HTML generation."""

    def test_report_contains_neighborhood_content(self, tmp_path):
        """Full report should include neighborhood-related HTML when
        ANI matrix and novel reads are present."""
        from metadarkmatter.visualization.report.generator import ReportGenerator

        ani, genomes = _make_multi_genus_ani_matrix()
        df = _make_novel_classifications(genomes)

        # Build a Polars ANI matrix suitable for ReportGenerator
        ani_rows: list[dict] = {"genome": genomes}
        for g1 in genomes:
            ani_rows[g1] = [ani.get_ani(g1, g2) for g2 in genomes]
        ani_pl = pl.DataFrame(ani_rows)

        generator = ReportGenerator(
            classifications=df,
            ani_matrix=ani_pl,
        )

        output_path = tmp_path / "report_neighborhood.html"
        generator.generate(output_path)

        assert output_path.exists()
        content = output_path.read_text()

        # The report should contain the novel diversity tab
        assert 'id="novel-diversity"' in content

        # Neighborhood-specific content should appear in the cluster table
        assert "Phylogenetic Context" in content or "neighborhood" in content.lower()

    def test_novel_section_cluster_table_has_neighborhood_columns(self):
        """The cluster table HTML built from enriched clusters should
        contain neighborhood context and placement support columns."""
        from metadarkmatter.visualization.report.novel_section import (
            build_cluster_table_html,
        )

        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        enriched = nbr_analyzer.analyze(clusters)

        html = build_cluster_table_html(enriched)

        # Table should contain neighborhood panel rows
        assert "neighborhood-panel" in html
        # Table should mention placement support
        assert "Placement" in html or "Support" in html
        # Each enriched cluster should have a genera sub-table
        assert "nbr-genera-table" in html

    def test_report_no_crash_without_novel_reads(self, tmp_path):
        """Report generation should complete without errors when
        there are no novel reads to cluster."""
        from metadarkmatter.visualization.report.generator import ReportGenerator

        _, genomes = _make_multi_genus_ani_matrix()

        # Only Known Species reads
        rows = []
        for i in range(50):
            rows.append({
                "read_id": f"known_{i}",
                "best_match_genome": genomes[0],
                "top_hit_identity": 98.0,
                "novelty_index": 2.0,
                "placement_uncertainty": 0.5,
                "taxonomic_call": "Known Species",
                "num_ambiguous_hits": 1,
            })
        df = pl.DataFrame(rows)

        generator = ReportGenerator(classifications=df)

        output_path = tmp_path / "report_known_only.html"
        generator.generate(output_path)

        assert output_path.exists()
        content = output_path.read_text()

        # Novel diversity tab should be absent
        assert 'id="novel-diversity"' not in content


class TestClusteringToNeighborhoodRoundTrip:
    """Verify data integrity across the clustering-neighborhood boundary."""

    def test_cluster_ids_preserved(self):
        """Cluster IDs should be identical before and after neighborhood
        enrichment."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()
        original_ids = {c.cluster_id for c in clusters}

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        enriched = nbr_analyzer.analyze(clusters)
        enriched_ids = {c.cluster_id for c in enriched}

        assert original_ids == enriched_ids

    def test_read_counts_preserved(self):
        """Read counts should not change through neighborhood analysis."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()
        original_counts = {c.cluster_id: c.read_count for c in clusters}

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        enriched = nbr_analyzer.analyze(clusters)
        enriched_counts = {c.cluster_id: c.read_count for c in enriched}

        assert original_counts == enriched_counts

    def test_neighborhood_cluster_id_matches_parent(self):
        """The cluster_id in each neighborhood should match its parent
        cluster's ID."""
        ani, genomes = _make_multi_genus_ani_matrix()
        genus_map = _make_genus_map(genomes)
        df = _make_novel_classifications(genomes)

        analyzer = NovelDiversityAnalyzer(
            classifications=df,
            metadata=None,
            ani_matrix=ani,
            min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()

        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        enriched = nbr_analyzer.analyze(clusters)

        for cluster in enriched:
            assert cluster.neighborhood is not None
            assert cluster.neighborhood.cluster_id == cluster.cluster_id
