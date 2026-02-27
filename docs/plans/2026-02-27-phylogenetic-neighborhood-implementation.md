# Phylogenetic Neighborhood Analysis Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add phylogenetic neighborhood analysis to novel diversity clusters, including adaptive genus boundary detection, placement support scores, and phylogenetic context text.

**Architecture:** Extend the existing novel diversity pipeline with a new neighborhood analyzer module. Genus boundary detection extends the adaptive.py GMM approach. Each novel cluster gets a PhylogeneticNeighborhood attached after clustering. Report enriched with new columns and collapsible panels.

**Tech Stack:** Python, Polars, Pydantic, scikit-learn (optional, for GMM), Plotly (for scatter boundary lines)

---

### Task 1: Add data models for PhylogeneticNeighborhood

**Files:**
- Modify: `src/metadarkmatter/core/novel_diversity/models.py`
- Test: `tests/unit/test_novel_clustering.py`

**Step 1: Write the failing test**

```python
# In tests/unit/test_novel_clustering.py
from metadarkmatter.core.novel_diversity.models import (
    GenusDistance,
    PhylogeneticNeighborhood,
    NovelCluster,
)

class TestPhylogeneticNeighborhood:
    def test_genus_distance_creation(self):
        gd = GenusDistance(
            genus="Francisella",
            representative_genome="GCF_000195955.2",
            estimated_ani=78.5,
            num_genomes_in_genus=5,
        )
        assert gd.genus == "Francisella"
        assert gd.estimated_ani == 78.5

    def test_neighborhood_creation(self):
        gd1 = GenusDistance(
            genus="Francisella",
            representative_genome="GCF_001",
            estimated_ani=78.5,
            num_genomes_in_genus=5,
        )
        gd2 = GenusDistance(
            genus="Allofrancisella",
            representative_genome="GCF_002",
            estimated_ani=72.0,
            num_genomes_in_genus=2,
        )
        neighborhood = PhylogeneticNeighborhood(
            cluster_id="NGN_001",
            nearest_genera=[gd1, gd2],
            placement_support=85.0,
            isolation_score=6.5,
            neighborhood_density=2,
            phylogenetic_context="Sister to Francisella (78% ANI), 6% isolated from Allofrancisella. Support: 85/100.",
            genus_boundary_ani=82.0,
        )
        assert neighborhood.placement_support == 85.0
        assert len(neighborhood.nearest_genera) == 2
        assert neighborhood.nearest_genera[0].genus == "Francisella"

    def test_novel_cluster_with_neighborhood(self):
        """NovelCluster should accept an optional neighborhood field."""
        cluster = NovelCluster(
            cluster_id="NGN_001",
            taxonomic_call="Novel Genus",
            nearest_genome="GCF_001",
            nearest_species="F. tularensis",
            nearest_genus="Francisella",
            nearest_family="Francisellaceae",
            novelty_band=20,
            read_count=47,
            mean_novelty_index=21.5,
            novelty_min=20.1,
            novelty_max=24.3,
            mean_placement_uncertainty=1.2,
            suggested_name="Francisellaceae gen. nov. MDM-001",
            confidence="High",
            phylogenetic_context="Novel genus within Francisellaceae",
        )
        # neighborhood is None by default
        assert cluster.neighborhood is None
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/test_novel_clustering.py::TestPhylogeneticNeighborhood -v`
Expected: FAIL (GenusDistance and PhylogeneticNeighborhood not defined, NovelCluster has no `neighborhood` field)

**Step 3: Write minimal implementation**

In `src/metadarkmatter/core/novel_diversity/models.py`, add:

```python
class GenusDistance(BaseModel):
    """Distance from a novel cluster to a reference genus."""
    genus: str = Field(description="Genus name")
    representative_genome: str = Field(description="Closest genome in this genus")
    estimated_ani: float = Field(description="Estimated ANI from cluster to this genus")
    num_genomes_in_genus: int = Field(description="Reference genomes in this genus")

    model_config = {"frozen": True}


class PhylogeneticNeighborhood(BaseModel):
    """Phylogenetic neighborhood profile for a novel cluster."""
    cluster_id: str = Field(description="Cluster this neighborhood belongs to")
    nearest_genera: list[GenusDistance] = Field(
        description="Nearest genera sorted by ANI (descending), top 5"
    )
    placement_support: float = Field(
        ge=0, le=100,
        description="Placement support score (0-100)"
    )
    isolation_score: float = Field(
        ge=0,
        description="ANI gap between nearest and second-nearest genus"
    )
    neighborhood_density: int = Field(
        ge=0,
        description="Number of genera within genus_boundary + 5% ANI"
    )
    phylogenetic_context: str = Field(
        description="One-line human-readable placement text"
    )
    genus_boundary_ani: float | None = Field(
        default=None,
        description="Detected genus boundary ANI from GMM"
    )

    model_config = {"frozen": True}
```

Add to `NovelCluster`:

```python
    neighborhood: PhylogeneticNeighborhood | None = Field(
        default=None,
        description="Phylogenetic neighborhood profile (computed post-clustering)"
    )
```

Note: `NovelCluster` has `model_config = {"frozen": True}`, so we need to add `neighborhood` as a regular field with a default of `None`. This means clusters are still created without neighborhood, and we'll construct new clusters with the neighborhood attached using `model_copy(update=...)` after neighborhood analysis.

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/test_novel_clustering.py::TestPhylogeneticNeighborhood -v`
Expected: PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/core/novel_diversity/models.py tests/unit/test_novel_clustering.py
git commit -m "feat: add GenusDistance and PhylogeneticNeighborhood data models"
```

---

### Task 2: Add genus boundary detection to adaptive.py

**Files:**
- Modify: `src/metadarkmatter/core/classification/adaptive.py`
- Test: `tests/unit/test_adaptive.py`

**Step 1: Write the failing test**

```python
# In tests/unit/test_adaptive.py
from metadarkmatter.core.classification.adaptive import (
    AdaptiveGenusThreshold,
    detect_genus_boundary,
)

class TestDetectGenusBoundary:
    def test_fallback_with_few_genomes(self):
        """Should fall back to default 80% when too few genomes."""
        # Create a small ANI matrix with < 5 genomes
        ani = _make_ani_matrix(["GCF_001", "GCF_002", "GCF_003"])
        result = detect_genus_boundary(ani)
        assert result.method == "fallback"
        assert result.genus_boundary == 80.0
        assert result.novelty_genus_min == 20.0

    def test_genus_boundary_reasonable_range(self):
        """Detected boundary should be within 70-90% ANI."""
        # Create an ANI matrix with clear genus-level structure:
        # 3 genera, 3 genomes each
        # Within-genus: 85-95% ANI, Between-genus: 70-80% ANI
        genomes = [f"G{i}" for i in range(9)]
        ani = _make_multi_genus_ani_matrix(genomes, genera_count=3)
        result = detect_genus_boundary(ani)
        assert 70.0 <= result.genus_boundary <= 90.0

    def test_genus_boundary_below_species_boundary(self):
        """Genus boundary must be lower than species boundary."""
        genomes = [f"G{i}" for i in range(9)]
        ani = _make_multi_genus_ani_matrix(genomes, genera_count=3)
        genus_result = detect_genus_boundary(ani)
        species_result = detect_species_boundary(ani)
        if genus_result.method != "fallback" and species_result.method != "fallback":
            assert genus_result.genus_boundary < species_result.species_boundary

    def test_metadata_validation(self):
        """When metadata provided, compute empirical inter-genus ANI range."""
        genomes = [f"G{i}" for i in range(9)]
        ani = _make_multi_genus_ani_matrix(genomes, genera_count=3)
        # Create mock metadata with genus labels
        genus_map = {f"G{i}": f"Genus_{i // 3}" for i in range(9)}
        result = detect_genus_boundary(ani, genus_map=genus_map)
        assert result.inter_genus_ani_range is not None
```

You'll need a helper `_make_multi_genus_ani_matrix` that creates an ANI matrix with clear within-genus (85-95%) and between-genus (70-80%) structure. Pattern it after the existing `_make_ani_matrix` helper in the test file.

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/test_adaptive.py::TestDetectGenusBoundary -v`
Expected: FAIL (AdaptiveGenusThreshold and detect_genus_boundary not defined)

**Step 3: Write minimal implementation**

In `src/metadarkmatter/core/classification/adaptive.py`, add:

```python
@dataclass
class AdaptiveGenusThreshold:
    """Result of adaptive genus boundary detection."""
    genus_boundary: float           # Detected ANI boundary (e.g., 82%)
    novelty_genus_min: float        # = 100 - genus_boundary
    confidence: float               # Separation quality 0-1
    method: str                     # "gmm_3component", "metadata_empirical", "fallback"
    inter_genus_ani_range: tuple[float, float] | None  # Observed min-max from metadata
    ani_values_used: int


def detect_genus_boundary(
    ani_matrix: ANIMatrix,
    min_genomes: int = 8,
    genus_map: dict[str, str] | None = None,
) -> AdaptiveGenusThreshold:
    """
    Detect natural genus boundary from ANI matrix distribution.

    Fits a 3-component GMM to pairwise ANI values. The three components
    typically correspond to:
    - Within-species comparisons (high ANI)
    - Between-species-within-genus comparisons (medium ANI)
    - Between-genera comparisons (low ANI)

    The genus boundary is estimated as the crossing point between the
    medium and low components.

    If genus_map is provided (genome -> genus name), compute empirical
    inter-genus ANI range for validation.

    Falls back to 80% ANI (20% novelty) if GMM fails.
    """
    default_boundary = 80.0
    default_novelty = 20.0
    genomes = sorted(ani_matrix.genomes)

    if len(genomes) < min_genomes:
        logger.info(
            "Too few genomes (%d < %d) for genus boundary detection. "
            "Using default boundary %.1f%%.",
            len(genomes), min_genomes, default_boundary,
        )
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=default_novelty,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=None,
            ani_values_used=0,
        )

    # Extract pairwise ANI values
    ani_values = []
    for i, g1 in enumerate(genomes):
        for j in range(i + 1, len(genomes)):
            g2 = genomes[j]
            ani = ani_matrix.get_ani(g1, g2)
            if ani > 0:
                ani_values.append(ani)

    # Compute empirical inter-genus range if genus labels available
    inter_genus_range = None
    if genus_map:
        inter_genus_anis = []
        for i, g1 in enumerate(genomes):
            for j in range(i + 1, len(genomes)):
                g2 = genomes[j]
                genus1 = genus_map.get(g1)
                genus2 = genus_map.get(g2)
                if genus1 and genus2 and genus1 != genus2:
                    ani = ani_matrix.get_ani(g1, g2)
                    if ani > 0:
                        inter_genus_anis.append(ani)
        if inter_genus_anis:
            inter_genus_range = (min(inter_genus_anis), max(inter_genus_anis))
            logger.info(
                "Empirical inter-genus ANI range: %.1f%% - %.1f%% (n=%d pairs).",
                inter_genus_range[0], inter_genus_range[1], len(inter_genus_anis),
            )

    if len(ani_values) < 10:
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=default_novelty,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=inter_genus_range,
            ani_values_used=len(ani_values),
        )

    ani_array = np.array(ani_values).reshape(-1, 1)

    try:
        from sklearn.mixture import GaussianMixture

        gmm = GaussianMixture(
            n_components=3,
            covariance_type="full",
            max_iter=300,
            n_init=10,
            random_state=42,
        )
        gmm.fit(ani_array)

        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())

        # Sort by mean: lowest ANI first
        order = np.argsort(means)
        sorted_means = means[order]
        sorted_stds = stds[order]

        # Genus boundary = crossing between component 0 (between-genera)
        # and component 1 (between-species-within-genus)
        low_mean, mid_mean = sorted_means[0], sorted_means[1]
        low_std, mid_std = sorted_stds[0], sorted_stds[1]

        # Separation quality between lowest two components
        separation = abs(mid_mean - low_mean) / (low_std + mid_std + 1e-10)

        if separation < 0.8:
            logger.info(
                "GMM genus components not well separated (%.2f). Using fallback.",
                separation,
            )
            return AdaptiveGenusThreshold(
                genus_boundary=default_boundary,
                novelty_genus_min=default_novelty,
                confidence=separation / 3.0,
                method="fallback",
                inter_genus_ani_range=inter_genus_range,
                ani_values_used=len(ani_values),
            )

        # Weighted midpoint between low and mid components
        w_low = 1.0 / (low_std + 1e-10)
        w_mid = 1.0 / (mid_std + 1e-10)
        boundary = (low_mean * w_low + mid_mean * w_mid) / (w_low + w_mid)

        # Clamp to reasonable genus range
        boundary = max(70.0, min(90.0, boundary))
        novelty_min = 100.0 - boundary
        confidence = min(1.0, separation / 3.0)

        # Validate against empirical range if available
        if inter_genus_range:
            empirical_mid = (inter_genus_range[0] + inter_genus_range[1]) / 2.0
            if abs(boundary - empirical_mid) > 5.0:
                logger.warning(
                    "GMM genus boundary (%.1f%%) diverges from empirical (%.1f%%). "
                    "Consider checking metadata genus labels.",
                    boundary, empirical_mid,
                )

        logger.info(
            "Detected genus boundary at %.1f%% ANI (novelty_genus_min=%.1f%%, "
            "confidence=%.2f, method=gmm_3component, n=%d).",
            boundary, novelty_min, confidence, len(ani_values),
        )

        return AdaptiveGenusThreshold(
            genus_boundary=boundary,
            novelty_genus_min=novelty_min,
            confidence=confidence,
            method="gmm_3component",
            inter_genus_ani_range=inter_genus_range,
            ani_values_used=len(ani_values),
        )

    except ImportError:
        logger.warning(
            "scikit-learn not available for genus boundary detection. "
            "Install with: pip install metadarkmatter[adaptive]"
        )
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=default_novelty,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=inter_genus_range,
            ani_values_used=len(ani_values),
        )
    except Exception as e:
        logger.warning("GMM genus boundary fitting failed: %s. Using fallback.", e)
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=default_novelty,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=inter_genus_range,
            ani_values_used=len(ani_values),
        )
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/test_adaptive.py::TestDetectGenusBoundary -v`
Expected: PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/core/classification/adaptive.py tests/unit/test_adaptive.py
git commit -m "feat: add genus boundary detection via 3-component GMM"
```

---

### Task 3: Create PhylogeneticNeighborhoodAnalyzer

**Files:**
- Create: `src/metadarkmatter/core/novel_diversity/neighborhood.py`
- Test: `tests/unit/test_neighborhood.py`

**Step 1: Write the failing test**

Create `tests/unit/test_neighborhood.py`:

```python
"""Tests for phylogenetic neighborhood analysis."""
from __future__ import annotations

import numpy as np
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.novel_diversity.models import NovelCluster
from metadarkmatter.core.novel_diversity.neighborhood import (
    PhylogeneticNeighborhoodAnalyzer,
)


def _make_ani_matrix_with_genera(
    genus_genomes: dict[str, list[str]],
    within_genus_ani: float = 90.0,
    between_genus_ani: float = 75.0,
) -> ANIMatrix:
    """Create ANI matrix with genus structure."""
    all_genomes = []
    genome_to_genus = {}
    for genus, genomes in genus_genomes.items():
        all_genomes.extend(genomes)
        for g in genomes:
            genome_to_genus[g] = genus

    n = len(all_genomes)
    arr = np.full((n, n), between_genus_ani, dtype=np.float32)
    np.fill_diagonal(arr, 100.0)

    for i in range(n):
        for j in range(n):
            if i != j and genome_to_genus[all_genomes[i]] == genome_to_genus[all_genomes[j]]:
                arr[i, j] = within_genus_ani

    matrix = ANIMatrix.__new__(ANIMatrix)
    matrix._genomes = tuple(all_genomes)
    matrix._genome_to_idx = {g: i for i, g in enumerate(all_genomes)}
    matrix._ani_array = arr
    matrix._default_ani = 70.0
    matrix._num_genomes = n
    return matrix


def _make_novel_cluster(
    cluster_id: str = "NGN_001",
    taxonomic_call: str = "Novel Genus",
    nearest_genome: str = "GCF_001",
    nearest_genus: str = "Francisella",
    mean_novelty_index: float = 21.5,
    read_count: int = 47,
    mean_bayesian_confidence: float = 75.0,
) -> NovelCluster:
    return NovelCluster(
        cluster_id=cluster_id,
        taxonomic_call=taxonomic_call,
        nearest_genome=nearest_genome,
        nearest_species="F. tularensis",
        nearest_genus=nearest_genus,
        nearest_family="Francisellaceae",
        novelty_band=20,
        read_count=read_count,
        mean_novelty_index=mean_novelty_index,
        novelty_min=20.1,
        novelty_max=24.3,
        mean_placement_uncertainty=1.2,
        suggested_name="Francisellaceae gen. nov. MDM-001",
        confidence="High",
        phylogenetic_context="Novel genus within Francisellaceae",
        mean_bayesian_confidence=mean_bayesian_confidence,
    )


class TestPhylogeneticNeighborhoodAnalyzer:
    def test_analyze_returns_neighborhoods(self):
        """Should return a neighborhood for each cluster."""
        genus_genomes = {
            "Francisella": ["GCF_001", "GCF_002"],
            "Allofrancisella": ["GCF_003"],
        }
        ani = _make_ani_matrix_with_genera(genus_genomes)
        genus_map = {}
        for genus, genomes in genus_genomes.items():
            for g in genomes:
                genus_map[g] = genus

        cluster = _make_novel_cluster(nearest_genome="GCF_001")
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani,
            genus_map=genus_map,
        )
        results = analyzer.analyze([cluster])

        assert len(results) == 1
        nbr = results[0].neighborhood
        assert nbr is not None
        assert nbr.cluster_id == "NGN_001"
        assert len(nbr.nearest_genera) >= 1
        assert nbr.nearest_genera[0].genus == "Francisella"  # Nearest
        assert 0 <= nbr.placement_support <= 100

    def test_isolation_score_computed(self):
        """Isolation score = gap between 1st and 2nd nearest genus."""
        genus_genomes = {
            "Francisella": ["GCF_001", "GCF_002"],
            "Allofrancisella": ["GCF_003"],
            "Thioalkalivibrio": ["GCF_004"],
        }
        ani = _make_ani_matrix_with_genera(
            genus_genomes, within_genus_ani=90.0, between_genus_ani=75.0,
        )
        genus_map = {}
        for genus, genomes in genus_genomes.items():
            for g in genomes:
                genus_map[g] = genus

        cluster = _make_novel_cluster(nearest_genome="GCF_001")
        analyzer = PhylogeneticNeighborhoodAnalyzer(ani_matrix=ani, genus_map=genus_map)
        results = analyzer.analyze([cluster])

        nbr = results[0].neighborhood
        assert nbr is not None
        # All between-genus ANI is 75%, so isolation between genera = 0
        assert nbr.isolation_score >= 0

    def test_phylogenetic_context_text_generated(self):
        """Should generate a human-readable context string."""
        genus_genomes = {
            "Francisella": ["GCF_001"],
            "Allofrancisella": ["GCF_002"],
        }
        ani = _make_ani_matrix_with_genera(genus_genomes)
        genus_map = {"GCF_001": "Francisella", "GCF_002": "Allofrancisella"}

        cluster = _make_novel_cluster(nearest_genome="GCF_001")
        analyzer = PhylogeneticNeighborhoodAnalyzer(ani_matrix=ani, genus_map=genus_map)
        results = analyzer.analyze([cluster])

        nbr = results[0].neighborhood
        assert nbr is not None
        assert "Francisella" in nbr.phylogenetic_context
        assert "ANI" in nbr.phylogenetic_context

    def test_novel_species_context(self):
        """Novel Species clusters should get species-level context."""
        genus_genomes = {
            "Francisella": ["GCF_001", "GCF_002"],
        }
        ani = _make_ani_matrix_with_genera(genus_genomes, within_genus_ani=92.0)
        genus_map = {"GCF_001": "Francisella", "GCF_002": "Francisella"}

        cluster = _make_novel_cluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_001",
            mean_novelty_index=8.0,
        )
        analyzer = PhylogeneticNeighborhoodAnalyzer(ani_matrix=ani, genus_map=genus_map)
        results = analyzer.analyze([cluster])

        nbr = results[0].neighborhood
        assert nbr is not None
        assert "Within" in nbr.phylogenetic_context or "within" in nbr.phylogenetic_context

    def test_without_metadata_uses_genome_groups(self):
        """Without genus labels, should use ANI-based genome neighborhoods."""
        genus_genomes = {
            "A": ["GCF_001", "GCF_002"],
            "B": ["GCF_003"],
        }
        ani = _make_ani_matrix_with_genera(genus_genomes)

        cluster = _make_novel_cluster(nearest_genome="GCF_001")
        # No genus_map
        analyzer = PhylogeneticNeighborhoodAnalyzer(ani_matrix=ani, genus_map=None)
        results = analyzer.analyze([cluster])

        nbr = results[0].neighborhood
        assert nbr is not None
        # Should still have nearest_genera, but with group names instead of genus names
        assert len(nbr.nearest_genera) >= 1
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/test_neighborhood.py -v`
Expected: FAIL (neighborhood.py doesn't exist)

**Step 3: Write implementation**

Create `src/metadarkmatter/core/novel_diversity/neighborhood.py`:

The analyzer should:
1. Accept ANI matrix, optional genus_map (genome->genus), optional genus_boundary
2. For each cluster:
   a. Use `cluster.nearest_genome` to estimate ANI to all reference genomes via triangulation (same formula as extended_matrix_builder: `100 - (divergence_novel + divergence_ref * 0.7)`)
   b. Group by genus (or by ANI neighborhoods via union-find at 80% if no genus_map)
   c. For each genus/group, take max ANI as the cluster-to-genus distance
   d. Sort by ANI descending, take top 5
   e. Compute isolation_score, neighborhood_density, placement_support
   f. Generate context text
3. Return clusters with neighborhood attached via `model_copy(update={"neighborhood": nbr})`

Full implementation is ~150-200 lines. Key methods:
- `analyze(clusters) -> list[NovelCluster]`
- `_compute_cluster_to_genome_anis(cluster) -> dict[str, float]`
- `_group_by_genus(genome_anis) -> list[GenusDistance]`
- `_compute_placement_support(cluster, nearest_genera, genus_boundary) -> float`
- `_generate_context_text(cluster, nearest_genera, placement_support) -> str`

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/test_neighborhood.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/core/novel_diversity/neighborhood.py tests/unit/test_neighborhood.py
git commit -m "feat: add PhylogeneticNeighborhoodAnalyzer for novel cluster context"
```

---

### Task 4: Update __init__.py exports

**Files:**
- Modify: `src/metadarkmatter/core/novel_diversity/__init__.py`

**Step 1: Add exports**

Add `PhylogeneticNeighborhoodAnalyzer`, `GenusDistance`, `PhylogeneticNeighborhood` to the module's `__init__.py`.

**Step 2: Verify import works**

Run: `python -c "from metadarkmatter.core.novel_diversity import PhylogeneticNeighborhoodAnalyzer; print('OK')"`

**Step 3: Commit**

```bash
git add src/metadarkmatter/core/novel_diversity/__init__.py
git commit -m "chore: export neighborhood analysis types from novel_diversity"
```

---

### Task 5: Integrate neighborhood analysis into report generator

**Files:**
- Modify: `src/metadarkmatter/visualization/report/generator.py` (the `_build_novel_diversity_section` method)

**Step 1: Write the failing test**

```python
# In tests/unit/test_report_generator.py
class TestNovelDiversityNeighborhood:
    def test_novel_section_includes_neighborhood_columns(self, tmp_path):
        """Cluster table should include Phylogenetic Context and Support columns."""
        # Create a report with novel reads and metadata
        # ... (use existing test fixtures pattern from the file)
        # Assert the HTML contains "Phylogenetic Context" and "Support" column headers
        pass
```

**Step 2: Modify `_build_novel_diversity_section`**

After clusters are formed (line ~1306), add:

```python
# Compute phylogenetic neighborhoods if ANI matrix available
if ani_matrix_obj is not None and self.genome_metadata is not None:
    from metadarkmatter.core.novel_diversity.neighborhood import (
        PhylogeneticNeighborhoodAnalyzer,
    )
    genus_map = {}
    for genome in ani_matrix_obj.genomes:
        genus = self.genome_metadata.get_genus(genome)
        if genus:
            genus_map[genome] = genus

    if genus_map:
        nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
            ani_matrix=ani_matrix_obj,
            genus_map=genus_map,
        )
        clusters = nbr_analyzer.analyze(clusters)
```

**Step 3: Run tests**

Run: `python -m pytest tests/unit/test_report_generator.py -v`

**Step 4: Commit**

```bash
git add src/metadarkmatter/visualization/report/generator.py tests/unit/test_report_generator.py
git commit -m "feat: integrate neighborhood analysis into report novel diversity section"
```

---

### Task 6: Enrich cluster table with neighborhood columns

**Files:**
- Modify: `src/metadarkmatter/visualization/report/novel_section.py`
- Modify: `src/metadarkmatter/visualization/report/styles.py`

**Step 1: Update cluster table template**

Add two new columns to `NOVEL_CLUSTER_TABLE_TEMPLATE`:
- "Context" column (after "Phylogenetic Placement")
- "Support" column (after "Confidence")

Update `NOVEL_CLUSTER_ROW_TEMPLATE` to include:
- `{neighborhood_context}` — from `cluster.neighborhood.phylogenetic_context` (or "-" if no neighborhood)
- `{placement_support}` — from `cluster.neighborhood.placement_support` (or "-")
- Support score with color-coded badge

**Step 2: Update the row building code**

In `novel_section.py`, the function that builds cluster rows needs to check `cluster.neighborhood` and extract values.

**Step 3: Add CSS for support score badges**

In `styles.py`:
```css
.support-high { color: #22c55e; }
.support-medium { color: #f59e0b; }
.support-low { color: #ef4444; }
```

**Step 4: Run tests**

Run: `python -m pytest tests/unit/test_report_generator.py -v`

**Step 5: Commit**

```bash
git add src/metadarkmatter/visualization/report/novel_section.py src/metadarkmatter/visualization/report/styles.py
git commit -m "feat: add neighborhood context and support columns to cluster table"
```

---

### Task 7: Add collapsible neighborhood panel per cluster

**Files:**
- Modify: `src/metadarkmatter/visualization/report/novel_section.py`
- Modify: `src/metadarkmatter/visualization/report/styles.py`

**Step 1: Create neighborhood panel template**

Add a `NEIGHBORHOOD_PANEL_TEMPLATE` that renders below each cluster row:
- Full list of nearest genera with ANI distances (mini table)
- Isolation score and component breakdown
- Placement support breakdown (each component with its contribution)
- Genus boundary and detection method

Use the existing `COLLAPSIBLE_PANEL_TEMPLATE` pattern from `templates.py`.

**Step 2: Wire into cluster table rendering**

After each cluster row, conditionally render the neighborhood panel if `cluster.neighborhood is not None`.

**Step 3: Add CSS for neighborhood panels**

Style the mini-table and score breakdown within the collapsible panel.

**Step 4: Run tests**

Run: `python -m pytest tests/unit/test_report_generator.py -v`

**Step 5: Commit**

```bash
git add src/metadarkmatter/visualization/report/novel_section.py src/metadarkmatter/visualization/report/styles.py
git commit -m "feat: add collapsible neighborhood panels for novel clusters"
```

---

### Task 8: Add genus boundary line to cluster scatter plot

**Files:**
- Modify: `src/metadarkmatter/visualization/report/novel_section.py` (the `build_cluster_scatter_figure` function)

**Step 1: Modify scatter builder to accept genus_boundary parameter**

The `build_cluster_scatter_figure` function should accept an optional `genus_boundary: float | None` parameter.

When provided, add a vertical dashed line at `100 - genus_boundary` on the novelty axis (x-axis), labeled "Genus boundary ({genus_boundary}% ANI)".

This sits alongside the existing species boundary line.

**Step 2: Pass genus_boundary from generator**

In `_build_novel_diversity_section`, if genus boundary was detected, pass it to the scatter builder.

**Step 3: Run tests**

Run: `python -m pytest tests/unit/test_report_generator.py -v`

**Step 4: Commit**

```bash
git add src/metadarkmatter/visualization/report/novel_section.py src/metadarkmatter/visualization/report/generator.py
git commit -m "feat: add adaptive genus boundary line to cluster scatter plot"
```

---

### Task 9: Add genus boundary to Summary tab Key Findings

**Files:**
- Modify: `src/metadarkmatter/visualization/report/generator.py` (the `_build_summary_section` method)

**Step 1: Store genus boundary result on ReportGenerator**

Add an optional `genus_boundary: AdaptiveGenusThreshold | None` attribute to the ReportGenerator (or pass it through config). Set it during the novel diversity section build.

**Step 2: Add Key Finding alert**

In `_build_summary_section`, if genus boundary was detected with method != "fallback", add an alert card:
`"Genus boundary detected at {X}% ANI ({method}, confidence: {Y})"`

**Step 3: Run tests**

Run: `python -m pytest tests/unit/test_report_generator.py -v`

**Step 4: Commit**

```bash
git add src/metadarkmatter/visualization/report/generator.py
git commit -m "feat: show adaptive genus boundary in Summary Key Findings"
```

---

### Task 10: Integration test — full pipeline with neighborhood

**Files:**
- Create or modify: `tests/integration/test_neighborhood_integration.py`

**Step 1: Write integration test**

Create a test that:
1. Creates an ANI matrix with 3 genera, 3 genomes each
2. Creates a classification DataFrame with novel reads
3. Creates metadata with genus labels
4. Runs the full pipeline: NovelDiversityAnalyzer → PhylogeneticNeighborhoodAnalyzer
5. Verifies clusters have neighborhoods attached
6. Generates a report and verifies HTML contains neighborhood content

**Step 2: Run test**

Run: `python -m pytest tests/integration/test_neighborhood_integration.py -v`

**Step 3: Run full test suite**

Run: `python -m pytest tests/ -x -q`
Expected: All tests pass

**Step 4: Commit**

```bash
git add tests/integration/test_neighborhood_integration.py
git commit -m "test: add integration test for phylogenetic neighborhood pipeline"
```

---

### Task 11: Update CLI to pass genus boundary

**Files:**
- Modify: `src/metadarkmatter/cli/report.py`

**Step 1: Wire genus boundary detection into report CLI**

In the `generate` command, after loading the ANI matrix and metadata:
- If both are available, call `detect_genus_boundary()` with a genus_map built from metadata
- Store the result for passing to ReportGenerator

This is only needed if the genus boundary should be computed at CLI level rather than inside the report generator. Check whether the generator already has access to the ANI matrix (it does — `self.ani_matrix`), so the computation can stay inside `_build_novel_diversity_section`.

If keeping it inside the generator: this task becomes a no-op (just verify the integration works end-to-end).

**Step 2: Run full test suite**

Run: `python -m pytest tests/ -x -q`

**Step 3: Commit (if changes made)**

```bash
git add src/metadarkmatter/cli/report.py
git commit -m "feat: wire genus boundary detection into report generation"
```
