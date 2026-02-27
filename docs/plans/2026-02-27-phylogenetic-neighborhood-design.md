# Phylogenetic Neighborhood Analysis for Novel Diversity

**Date:** 2026-02-27
**Status:** Approved

## Problem

Novel Genus and Novel Species clusters are currently assigned through threshold arithmetic alone (20-25% novelty for genus, 4-20% for species). There is no consideration of the phylogenetic neighborhood — how many genera are nearby, whether the cluster sits on an isolated branch or in a densely sampled region, or what the actual genus boundary is for this particular family.

This leads to:
- Fixed thresholds that don't adapt to lineage-specific evolutionary rates
- Generic names ("Francisellaceae gen. nov. MDM-001") with no phylogenetic context
- Single-nearest-genome placement that ignores the broader neighborhood
- No placement confidence metric to help users assess reliability

## Design Goals

- Detect family-specific genus boundaries from the ANI distribution
- Compute a phylogenetic neighborhood profile for each novel cluster
- Generate placement support scores combining isolation, boundary, and evidence
- Produce phylogenetic context text for human interpretation
- Add collapsible neighborhood panels to the report for detailed exploration
- Apply to both Novel Species and Novel Genus clusters

## Architecture

### Component 1: Genus Boundary Detection

**File:** `src/metadarkmatter/core/classification/adaptive.py`

Extends the existing species boundary detection with a genus-level boundary.

**Algorithm:**
1. Extract all pairwise ANI values from the upper triangle of the ANI matrix
2. Fit a 3-component GMM: within-species, between-species-within-genus, between-genera
3. Genus boundary = crossing point between components 2 and 3
4. Validation: genus boundary must be lower than species boundary and within [70%, 90%]
5. If metadata genus labels are available, compute empirical inter-genus ANI range as validation
6. Fallback: use default 80% ANI if GMM fails and no metadata available
7. Log warning when metadata is absent (metadata strongly recommended)

**Data model:**
```python
@dataclass
class AdaptiveGenusThreshold:
    genus_boundary: float                          # Detected ANI boundary (e.g., 82%)
    novelty_genus_min: float                       # = 100 - genus_boundary (e.g., 18%)
    confidence: float                              # Separation quality 0-1
    method: str                                    # "gmm_3component", "metadata_empirical", "fallback"
    inter_genus_ani_range: tuple[float, float] | None  # Observed min-max from metadata
```

### Component 2: Phylogenetic Neighborhood Analyzer

**File:** `src/metadarkmatter/core/novel_diversity/neighborhood.py` (new)

Computes a neighborhood profile for each novel cluster.

**Data models** (in `core/novel_diversity/models.py`):
```python
@dataclass
class GenusDistance:
    genus: str                    # e.g., "Francisella"
    representative_genome: str    # Closest genome in that genus
    estimated_ani: float          # ANI from cluster to that genome
    num_genomes_in_genus: int     # How many genomes of this genus in the reference set

@dataclass
class PhylogeneticNeighborhood:
    cluster_id: str
    nearest_genera: list[GenusDistance]       # Sorted by ANI (descending), top 5
    placement_support: float                 # 0-100 score
    isolation_score: float                   # ANI gap between 1st and 2nd nearest genus
    neighborhood_density: int                # Genera within genus_boundary + 5% ANI
    phylogenetic_context: str                # One-line human-readable placement text
    genus_boundary_ani: float | None         # Detected genus boundary from GMM
```

**Algorithm per cluster:**
1. Compute ANI to all reference genomes using triangulation from the cluster's nearest genome
2. Group references by genus (from metadata). For each genus, take max ANI across its genomes
3. Compute isolation score: `ANI(nearest_genus) - ANI(second_nearest_genus)`
4. Compute neighborhood density: count genera within `genus_boundary + 5%` ANI
5. Compute placement support score (see formula below)
6. Generate one-line phylogenetic context text

**Without metadata:** Group genomes by ANI neighborhoods (union-find at 80%) and report genome groups instead of named genera.

### Placement Support Score

```
isolation_component   = min(isolation_score / 10.0, 1.0) * 40
boundary_component    = min((genus_boundary - estimated_ani_to_nearest) / 10.0, 1.0) * 30
read_component        = min(cluster.read_count / 20, 1.0) * 15
confidence_component  = (cluster.mean_bayesian_confidence / 100.0) * 15
placement_support     = isolation + boundary + read + confidence
```

- **Isolation (40%):** Gap between nearest and second-nearest genus. 10%+ gap = max score.
- **Boundary clarity (30%):** How far below the detected genus boundary. 10%+ below = max score.
- **Read support (15%):** Cluster size. 20+ reads = max score.
- **Bayesian confidence (15%):** From the classifier's entropy-based confidence.

### Phylogenetic Context Text

One-line format for table cells:

- Novel Genus, high support: `"Sister to Francisella (78% ANI), 6% isolated from Allofrancisella. Support: 85/100."`
- Novel Genus, low support: `"Between Francisella (79% ANI) and Allofrancisella (78% ANI). Support: 35/100."`
- Novel Species: `"Within Francisella, closest to F. tularensis (92% ANI), distinct from F. novicida (88% ANI)."`

### Component 3: Report Enrichment

**Cluster table** — Two new columns:
- **Phylogenetic Context**: One-line context text
- **Support**: 0-100 with color coding (green >= 70, orange 40-70, red < 40)

**Collapsible neighborhood panel** (per cluster):
- Full list of nearest genera with ANI distances
- Isolation score breakdown
- Placement support component breakdown
- Genus boundary used and detection method

**Cluster scatter plot:**
- Add detected genus boundary as a vertical dashed line
- Existing species boundary line remains

**Summary tab Key Findings:**
- If genus boundary adaptively detected: `"Genus boundary detected at 82% ANI (GMM, confidence: 0.85)"`

**Novel Diversity KPI strip:**
- Add genus boundary metric card

### Component 4: Pipeline Integration

1. `cli/score.py` — Call `detect_genus_boundary()` alongside `detect_species_boundary()`. Pass genus boundary to clustering step.

2. `core/novel_diversity/clustering.py` — After clusters formed, call `PhylogeneticNeighborhoodAnalyzer.analyze()`. Attach `PhylogeneticNeighborhood` to each cluster.

3. `core/novel_diversity/models.py` — Extend `NovelCluster` with `neighborhood: PhylogeneticNeighborhood | None`.

4. `cli/report.py` — Pass enriched cluster list to `ReportGenerator`.

5. `visualization/report/novel_section.py` — Render enriched cluster table + neighborhood panels.

**Backward compatibility:** If `neighborhood` is None (no metadata, old data), the report displays the current content unchanged.

## Files to Create/Modify

| File | Action |
|------|--------|
| `core/classification/adaptive.py` | Add `detect_genus_boundary()` |
| `core/novel_diversity/models.py` | Add `GenusDistance`, `PhylogeneticNeighborhood`, extend `NovelCluster` |
| `core/novel_diversity/neighborhood.py` | **New** — `PhylogeneticNeighborhoodAnalyzer` |
| `core/novel_diversity/clustering.py` | Integrate neighborhood analysis after clustering |
| `cli/score.py` | Call genus boundary detection, pass to clustering |
| `cli/report.py` | Pass enriched clusters to report |
| `visualization/report/novel_section.py` | Enriched cluster table + neighborhood panels |
| `visualization/report/styles.py` | Styles for support score badges and neighborhood panels |
| Tests | Unit tests for GMM genus detection, neighborhood analysis, placement support |

## What Stays Unchanged

- Bayesian 2D Gaussian classification (novelty x uncertainty)
- Stage 2 refinement logic
- Novel cluster ID naming scheme (NGN_###, NSP_###)
- Extended ANI heatmap builder
- Tree placement algorithm (branch lengths)
- Data tab and its filter chips
