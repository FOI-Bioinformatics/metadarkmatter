# Phylogenetic Tree Visualization Design

**Date:** 2026-02-05
**Status:** Draft
**Author:** Brainstorming session

## Overview

Add an interactive phylogenetic tree visualization to the metadarkmatter HTML report. The tree shows reference genomes and novel clusters in evolutionary context, enabling researchers to understand where novel diversity fits relative to known taxa.

## Goals

1. **Answer the core question:** "Where does this novel thing fit evolutionarily?"
2. **Multi-scale navigation:** Start at family level, drill into genera, then species neighborhoods
3. **Rich annotations:** Show novelty metrics, read counts, and confidence directly on tree nodes
4. **Flexible input:** Accept user-provided Newick trees or auto-generate from ANI matrix
5. **Non-disruptive:** Add as new tab; existing visualizations remain unchanged

## Non-Goals

- Replacing existing heatmaps or scatter plots
- Building a full phylogenetic inference pipeline (use existing tools for that)
- Cross-linking tree selections with other tabs (future enhancement)

---

## Architecture

### New Components

```
src/metadarkmatter/
├── core/
│   └── phylogeny/
│       ├── __init__.py
│       ├── tree_builder.py      # ANI -> neighbor-joining tree
│       └── placement.py         # Insert novel clusters into tree
├── visualization/
│   └── report/
│       ├── generator.py         # Add _build_phylogeny_section()
│       └── assets/
│           └── phylotree_bundle.js  # Phylotree.js + custom styling
```

### Data Flow

```
ANI Matrix ─────┬──> tree_builder.py ──> Newick string
                │                              │
User Newick ────┘                              v
                                        placement.py
                                              │
Novel Clusters ───────────────────────────────┘
                                              │
                                              v
                                    Tree with novel nodes
                                              │
                                              v
                                    _build_phylogeny_section()
                                              │
                                              v
                                    HTML report with Phylogeny tab
```

---

## Component Details

### 1. Tree Builder (`core/phylogeny/tree_builder.py`)

Converts ANI similarity matrix to a neighbor-joining phylogenetic tree.

```python
def ani_to_newick(ani_matrix: pd.DataFrame, genome_metadata: dict = None) -> str:
    """Convert ANI matrix to neighbor-joining tree in Newick format.

    Args:
        ani_matrix: Square DataFrame with ANI values (0-100 scale)
        genome_metadata: Optional dict mapping genome IDs to display names

    Returns:
        Newick-formatted tree string with branch lengths
    """
    # 1. Convert ANI (similarity) to distance
    distance_matrix = 100 - ani_matrix
    distance_matrix = distance_matrix.clip(lower=0, upper=50)  # Cap extreme distances

    # 2. Handle missing values (too divergent for ANI calculation)
    distance_matrix = distance_matrix.fillna(50.0)

    # 3. Build neighbor-joining tree using BioPython
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

    names = list(ani_matrix.columns)
    matrix = _to_lower_triangular(distance_matrix)
    dm = DistanceMatrix(names, matrix)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # 4. Root at midpoint (standard for unrooted trees)
    tree.root_at_midpoint()

    # 5. Export as Newick
    return tree.format("newick")


def load_user_tree(newick_path: Path, ani_matrix: pd.DataFrame) -> str:
    """Load and validate user-provided Newick tree.

    Args:
        newick_path: Path to Newick file
        ani_matrix: ANI matrix for validation

    Returns:
        Validated Newick string (possibly pruned)
    """
    tree = Phylo.read(newick_path, "newick")

    # Validate tip names against ANI matrix
    tree_tips = {tip.name for tip in tree.get_terminals()}
    matrix_genomes = set(ani_matrix.columns)

    missing = matrix_genomes - tree_tips
    extra = tree_tips - matrix_genomes

    if missing:
        logger.warning(f"{len(missing)} genomes in ANI matrix not in tree - excluded from phylogeny")
    if extra:
        tree = _prune_tips(tree, keep=matrix_genomes)
        logger.info(f"Pruned {len(extra)} unused tips from provided tree")

    return tree.format("newick")
```

**Key decisions:**

- **Distance cap at 50:** ANI below 50% is unreliable and would distort tree topology
- **Midpoint rooting:** Standard approach when no outgroup is specified
- **Graceful mismatch handling:** Warn and continue rather than fail

### 2. Novel Cluster Placement (`core/phylogeny/placement.py`)

Inserts novel clusters as leaf nodes at estimated positions.

```python
@dataclass
class NovelCluster:
    cluster_id: str
    classification: str  # "novel_species" or "novel_genus"
    best_match_genome: str
    mean_identity: float
    mean_novelty: float
    mean_uncertainty: float
    read_count: int
    confidence_rating: str  # "High", "Medium", "Low"


def place_novel_clusters(
    tree: Phylo.BaseTree,
    novel_clusters: list[NovelCluster],
    ani_matrix: pd.DataFrame,
    genome_metadata: pd.DataFrame = None
) -> Phylo.BaseTree:
    """Insert novel clusters as leaf nodes at estimated positions.

    Placement strategy:
    - Novel species: Insert as sibling to nearest reference genome
    - Novel genus: Insert as sibling to nearest genus clade

    Args:
        tree: Reference genome tree
        novel_clusters: List of novel clusters to place
        ani_matrix: For estimating distances
        genome_metadata: For genus-level placement decisions

    Returns:
        Tree with novel clusters inserted
    """
    for cluster in novel_clusters:
        # Estimate branch length from alignment identity
        estimated_distance = 100 - cluster.mean_identity

        if cluster.classification == "novel_genus" and genome_metadata is not None:
            # Place as sibling to entire genus clade
            target_clade = _find_genus_clade(tree, cluster.best_match_genome, genome_metadata)
        else:
            # Place as sibling to nearest reference
            target_clade = _find_clade_by_name(tree, cluster.best_match_genome)

        # Create novel node with metadata
        novel_node = Phylo.BaseTree.Clade(
            name=f"[NOVEL] {cluster.cluster_id}",
            branch_length=estimated_distance / 2
        )

        # Attach metadata for visualization
        novel_node.comment = json.dumps({
            "type": cluster.classification,
            "read_count": cluster.read_count,
            "mean_novelty": cluster.mean_novelty,
            "mean_uncertainty": cluster.mean_uncertainty,
            "confidence": cluster.confidence_rating,
            "nearest_ref": cluster.best_match_genome,
            "est_ani": cluster.mean_identity,
            "is_novel": True
        })

        # Insert as sibling
        _insert_as_sibling(tree, target_clade, novel_node, estimated_distance)

    return tree
```

**Insertion diagram:**

```
Before:                    After inserting Novel_X near Genome_A:

    +-- Genome_A               +-- Genome_A
----+                      ----+-- [new internal node]
    +-- Genome_B               |   +-- [NOVEL] Novel_X
                               +-- Genome_B
```

### 3. Visualization (`visualization/report/generator.py`)

New method to build the Phylogeny tab.

```python
def _build_phylogeny_section(
    self,
    classifications_df: pd.DataFrame,
    ani_matrix: pd.DataFrame,
    user_tree: Optional[Path] = None
) -> str:
    """Build interactive phylogenetic tree section.

    Returns:
        HTML string for Phylogeny tab content
    """
    # 1. Build or load tree
    if user_tree:
        newick = load_user_tree(user_tree, ani_matrix)
        tree_source = "user-provided"
    else:
        newick = ani_to_newick(ani_matrix)
        tree_source = "neighbor-joining from ANI"

    # 2. Extract novel clusters from classifications
    novel_clusters = self._extract_novel_clusters(classifications_df)

    # 3. Place novel clusters on tree
    tree = Phylo.read(StringIO(newick), "newick")
    if novel_clusters:
        tree = place_novel_clusters(tree, novel_clusters, ani_matrix, self.metadata)

    # 4. Prepare node annotations
    annotations = self._build_tree_annotations(tree, classifications_df)

    # 5. Determine initial collapse level for large trees
    genome_count = len(ani_matrix)
    initial_collapse = "genus" if genome_count > 100 else "none"

    # 6. Render HTML with embedded Phylotree.js
    return self._render_phylogeny_html(
        newick=tree.format("newick"),
        annotations=annotations,
        tree_source=tree_source,
        initial_collapse=initial_collapse,
        novel_count=len(novel_clusters)
    )
```

### 4. JavaScript Visualization

Using Phylotree.js for interactive tree rendering.

```javascript
// Embedded in report HTML template

document.addEventListener('DOMContentLoaded', function() {
    const container = document.getElementById('phylogeny-container');
    const newickString = TREE_DATA.newick;  // Injected by Python
    const annotations = TREE_DATA.annotations;

    const tree = new Phylotree(newickString);

    tree.render({
        container: container,
        width: 900,
        height: Math.max(600, TREE_DATA.tip_count * 18),

        // Layout
        layout: "rectangular",
        alignTips: true,

        // Interactivity
        collapsible: true,
        selectable: true,
        zoomable: true,

        // Custom styling
        nodeStyler: function(element, node) {
            const meta = annotations[node.name];

            if (meta && meta.is_novel) {
                // Novel clusters: colored circles sized by read count
                const color = meta.type === "novel_species" ? "#f59e0b" : "#ef4444";
                const radius = Math.sqrt(meta.read_count / 10) + 4;

                element.select("circle")
                    .attr("r", radius)
                    .attr("fill", color)
                    .attr("stroke", "#ffffff")
                    .attr("stroke-width", 2);

                // Add confidence badge
                element.append("text")
                    .attr("class", "confidence-badge")
                    .attr("dx", radius + 4)
                    .attr("dy", 4)
                    .text(meta.confidence[0])  // "H", "M", or "L"
                    .attr("fill", meta.confidence === "High" ? "#16a34a" :
                                  meta.confidence === "Medium" ? "#f59e0b" : "#6b7280");
            } else {
                // Reference genomes: small neutral circles
                element.select("circle")
                    .attr("r", 3)
                    .attr("fill", "#6b7280");
            }
        },

        branchStyler: function(element, node) {
            // Highlight branches leading to novel clusters
            if (node.name && node.name.startsWith("[NOVEL]")) {
                element.attr("stroke", "#f59e0b")
                       .attr("stroke-width", 2)
                       .attr("stroke-dasharray", "4,2");
            }
        }
    });

    // Tooltip on hover
    tree.on("node_hover", function(node, element, event) {
        const meta = annotations[node.name];
        if (meta && meta.is_novel) {
            showTooltip(event, `
                <strong>${node.name}</strong><br>
                Classification: ${meta.type.replace("_", " ")}<br>
                Reads: ${meta.read_count.toLocaleString()}<br>
                Mean Novelty: ${meta.mean_novelty.toFixed(1)}%<br>
                Mean Uncertainty: ${meta.mean_uncertainty.toFixed(1)}%<br>
                Confidence: ${meta.confidence}<br>
                Nearest Reference: ${meta.nearest_ref}<br>
                Est. ANI: ${meta.est_ani.toFixed(1)}%
            `);
        }
    });

    // Layout toggle
    document.getElementById('layout-toggle').addEventListener('click', function() {
        const currentLayout = tree.getLayout();
        tree.setLayout(currentLayout === "rectangular" ? "radial" : "rectangular");
        tree.render();
    });

    // Collapse controls
    document.getElementById('collapse-genera').addEventListener('click', function() {
        tree.collapseToLevel("genus");
    });

    document.getElementById('expand-all').addEventListener('click', function() {
        tree.expandAll();
    });

    // Initial collapse for large trees
    if (TREE_DATA.initial_collapse === "genus") {
        tree.collapseToLevel("genus");
    }
});
```

---

## CLI Interface

### New Options

```bash
# Auto-generate tree from ANI matrix (default)
metadarkmatter report generate \
    --classifications classifications.csv \
    --ani ani_matrix.csv \
    --output report.html

# Use provided phylogenetic tree
metadarkmatter report generate \
    --classifications classifications.csv \
    --ani ani_matrix.csv \
    --tree reference_tree.nwk \
    --output report.html

# Skip phylogeny tab (faster, smaller report)
metadarkmatter report generate \
    --classifications classifications.csv \
    --ani ani_matrix.csv \
    --no-phylogeny \
    --output report.html
```

### Parameter Definitions

```python
@app.command()
def generate(
    # ... existing parameters ...

    tree: Optional[Path] = typer.Option(
        None,
        "--tree", "-t",
        help="Newick tree file for phylogenetic context. "
             "If not provided, neighbor-joining tree is built from ANI matrix.",
        exists=True,
        dir_okay=False,
    ),
    no_phylogeny: bool = typer.Option(
        False,
        "--no-phylogeny",
        help="Skip phylogeny tab generation (faster for large datasets).",
    ),
):
```

---

## Edge Cases

### Small Datasets (< 3 genomes)

Neighbor-joining requires at least 3 taxa:

```python
if len(ani_matrix) < 3:
    logger.warning(
        "Too few genomes for phylogenetic tree (need >= 3). "
        "Skipping phylogeny tab."
    )
    return None
```

### Large Datasets (> 500 genomes)

Tree becomes unwieldy; start collapsed:

```python
if len(ani_matrix) > 500:
    logger.info("Large dataset: tree will start collapsed to genus level")
    initial_collapse = "genus"
```

### Missing ANI Values

Fill with maximum distance and warn:

```python
missing_count = ani_matrix.isna().sum().sum() // 2
if missing_count > 0:
    logger.warning(
        f"{missing_count} genome pairs lack ANI values; "
        "using maximum distance (50)"
    )
    distance_matrix = distance_matrix.fillna(50.0)
```

### Tree/Matrix Mismatch

When user-provided tree doesn't match ANI matrix:

- **Genomes in matrix but not tree:** Warn, exclude from phylogeny visualization
- **Tree tips not in matrix:** Prune from tree, log count

### No Novel Clusters

Tree still renders showing reference genome relationships:

```python
if not novel_clusters:
    logger.info("No novel clusters to place; showing reference tree only")
```

---

## Dependencies

### Python

No new dependencies required:

- `biopython` - Already present for sequence handling; provides `Bio.Phylo`
- `scipy` - Already present for heatmap clustering; provides distance calculations

### JavaScript

Phylotree.js bundled as static asset:

- **Source:** https://github.com/veg/phylotree.js
- **Size:** ~180KB minified
- **License:** MIT

Total report size increase: ~200-250KB (current reports are 2-5MB).

---

## Testing Strategy

### Unit Tests

```python
# tests/unit/core/phylogeny/test_tree_builder.py

def test_ani_to_newick_basic():
    """3x3 ANI matrix produces valid Newick."""
    ani = pd.DataFrame({
        "A": [100, 95, 80],
        "B": [95, 100, 82],
        "C": [80, 82, 100]
    }, index=["A", "B", "C"])

    newick = ani_to_newick(ani)

    assert newick.endswith(";")
    assert "A" in newick and "B" in newick and "C" in newick


def test_ani_to_newick_missing_values():
    """Missing ANI values handled gracefully."""
    ani = pd.DataFrame({
        "A": [100, 95, np.nan],
        "B": [95, 100, 82],
        "C": [np.nan, 82, 100]
    }, index=["A", "B", "C"])

    newick = ani_to_newick(ani)
    assert newick is not None


def test_ani_to_newick_minimum_size():
    """Fewer than 3 genomes returns None."""
    ani = pd.DataFrame({
        "A": [100, 95],
        "B": [95, 100]
    }, index=["A", "B"])

    result = ani_to_newick(ani)
    assert result is None


# tests/unit/core/phylogeny/test_placement.py

def test_place_novel_cluster_as_sibling():
    """Novel cluster inserted as sibling to nearest reference."""
    tree = Phylo.read(StringIO("((A:5,B:5):3,C:8);"), "newick")
    cluster = NovelCluster(
        cluster_id="Novel_1",
        classification="novel_species",
        best_match_genome="A",
        mean_identity=88.0,
        mean_novelty=12.0,
        mean_uncertainty=1.5,
        read_count=150,
        confidence_rating="High"
    )

    result = place_novel_clusters(tree, [cluster], ani_matrix)

    tip_names = [t.name for t in result.get_terminals()]
    assert "[NOVEL] Novel_1" in tip_names


def test_place_novel_genus_at_clade():
    """Novel genus placed as sibling to genus clade, not single species."""
    # ... test implementation
```

### Integration Tests

```python
# tests/integration/test_report_phylogeny.py

def test_report_generates_with_phylogeny_tab(tmp_path, sample_classifications):
    """Full report includes interactive phylogeny."""
    result = runner.invoke(app, [
        "report", "generate",
        "--classifications", str(sample_classifications.csv),
        "--ani", str(sample_classifications.ani),
        "--output", str(tmp_path / "report.html")
    ])

    assert result.exit_code == 0
    html = (tmp_path / "report.html").read_text()
    assert "Phylogeny" in html
    assert "phylotree" in html.lower()
    assert "newickString" in html or "TREE_DATA" in html


def test_report_with_user_tree(tmp_path, sample_classifications, sample_tree):
    """Report uses user-provided tree when specified."""
    result = runner.invoke(app, [
        "report", "generate",
        "--classifications", str(sample_classifications.csv),
        "--ani", str(sample_classifications.ani),
        "--tree", str(sample_tree),
        "--output", str(tmp_path / "report.html")
    ])

    assert result.exit_code == 0
    assert "user-provided" in result.output or result.exit_code == 0


def test_report_no_phylogeny_flag(tmp_path, sample_classifications):
    """--no-phylogeny skips phylogeny tab."""
    result = runner.invoke(app, [
        "report", "generate",
        "--classifications", str(sample_classifications.csv),
        "--ani", str(sample_classifications.ani),
        "--no-phylogeny",
        "--output", str(tmp_path / "report.html")
    ])

    assert result.exit_code == 0
    html = (tmp_path / "report.html").read_text()
    assert "Phylogeny" not in html
```

### Visual Regression Checklist

- [ ] Tree renders with correct topology for known test case
- [ ] Novel clusters appear with orange (species) or red (genus) markers
- [ ] Marker size scales with read count
- [ ] Confidence badges display correctly (H/M/L)
- [ ] Clade collapse/expand works via click
- [ ] "Collapse to Genera" and "Expand All" buttons function
- [ ] Hover tooltips show all metrics
- [ ] Zoom and pan work smoothly
- [ ] Radial layout toggle works
- [ ] Large trees (100+ tips) start collapsed appropriately

---

## Implementation Phases

### Phase 1: Core Tree Building
- Implement `tree_builder.py` with ANI-to-Newick conversion
- Implement `placement.py` for novel cluster insertion
- Unit tests for both modules

### Phase 2: Report Integration
- Add `_build_phylogeny_section()` to report generator
- Bundle Phylotree.js as static asset
- Basic tree rendering without annotations

### Phase 3: Rich Annotations
- Node styling by classification type
- Confidence badges
- Hover tooltips with full metrics
- Branch styling for novel clusters

### Phase 4: Interactivity
- Clade collapse/expand
- Layout toggle (rectangular/radial)
- Zoom and pan
- Collapse controls for large trees

### Phase 5: CLI and Polish
- Add `--tree` and `--no-phylogeny` CLI options
- Error handling for edge cases
- Integration tests
- Documentation updates

---

## Future Enhancements (Out of Scope)

- Cross-linking: Click tree node, highlight in scatter plots
- Export tree as SVG/PDF for publication
- Comparative view: Multiple samples on same tree
- Animation: Show how classifications change with different thresholds
