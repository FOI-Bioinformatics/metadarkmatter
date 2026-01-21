"""
Novel diversity section builder for HTML reports.

Creates visualizations and content for the Novel Diversity tab,
including cluster summary tables, phylogenetic context charts,
and cluster quality scatter plots.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import plotly.graph_objects as go
    from metadarkmatter.core.novel_diversity import (
        NovelCluster,
        NovelDiversitySummary,
    )

logger = logging.getLogger(__name__)


# =============================================================================
# HTML Templates for Novel Diversity Section
# =============================================================================

NOVEL_SUMMARY_TEMPLATE: str = '''
<div class="novel-summary">
    <div class="novel-intro">
        <h3>Novel Diversity Analysis</h3>
        <p>Reads classified as Novel Species or Novel Genus are grouped into clusters
           based on their nearest reference genome and divergence level. Each cluster
           represents a putative novel organism that may warrant further investigation.</p>
    </div>
    <div class="novel-metrics">
        <div class="novel-metric highlight">
            <div class="novel-metric-value">{total_clusters}</div>
            <div class="novel-metric-label">Novel Clusters</div>
            <div class="novel-metric-hint">{total_reads:,} reads total</div>
        </div>
        <div class="novel-metric species">
            <div class="novel-metric-value">{novel_species_clusters}</div>
            <div class="novel-metric-label">Novel Species</div>
            <div class="novel-metric-hint">{novel_species_reads:,} reads</div>
        </div>
        <div class="novel-metric genus">
            <div class="novel-metric-value">{novel_genus_clusters}</div>
            <div class="novel-metric-label">Novel Genus</div>
            <div class="novel-metric-hint">{novel_genus_reads:,} reads</div>
        </div>
        <div class="novel-metric confidence">
            <div class="novel-metric-value">{high_confidence}</div>
            <div class="novel-metric-label">High Confidence</div>
            <div class="novel-metric-hint">Priority for validation</div>
        </div>
    </div>
</div>
'''

NOVEL_CONFIDENCE_GUIDE_TEMPLATE: str = '''
<div class="novel-confidence-guide">
    <h4>Confidence Ratings</h4>
    <div class="confidence-grid">
        <div class="confidence-item high">
            <span class="conf-badge">High</span>
            <span class="conf-criteria">>=10 reads, uncertainty <5%, discovery >=75</span>
            <span class="conf-action">Prioritize for experimental validation</span>
        </div>
        <div class="confidence-item medium">
            <span class="conf-badge">Medium</span>
            <span class="conf-criteria">>=5 reads, uncertainty <10%, discovery >=50</span>
            <span class="conf-action">Include in candidate list</span>
        </div>
        <div class="confidence-item low">
            <span class="conf-badge">Low</span>
            <span class="conf-criteria">Other clusters meeting minimum size</span>
            <span class="conf-action">May need additional evidence</span>
        </div>
    </div>
</div>
'''

NOVEL_CLUSTER_TABLE_TEMPLATE: str = '''
<div class="novel-cluster-table">
    <h4>Novel Clusters</h4>
    <div class="table-wrapper">
        <table class="cluster-table">
            <thead>
                <tr>
                    <th>Cluster ID</th>
                    <th>Type</th>
                    <th>Reads</th>
                    <th>Phylogenetic Placement</th>
                    <th>Est. ANI</th>
                    <th>Uncertainty</th>
                    <th>Confidence</th>
                </tr>
            </thead>
            <tbody>
{cluster_rows}
            </tbody>
        </table>
    </div>
</div>
'''

NOVEL_CLUSTER_ROW_TEMPLATE: str = '''<tr class="cluster-row {confidence_class}">
    <td class="cluster-id">{cluster_id}</td>
    <td class="cluster-type {type_class}">{taxonomic_call}</td>
    <td class="cluster-reads">{read_count:,}</td>
    <td class="cluster-placement">
        <div class="placement-name">{suggested_name}</div>
        <div class="placement-context">{phylogenetic_context}</div>
        <div class="placement-ref" title="{full_taxonomy}">Nearest: {closest_taxon}</div>
    </td>
    <td class="cluster-ani">{estimated_ani:.1f}%</td>
    <td class="cluster-uncertainty">{mean_uncertainty:.1f}%</td>
    <td class="cluster-confidence"><span class="conf-tag {confidence_class}">{confidence}</span></td>
</tr>'''

NOVEL_PHYLOGENETIC_CONTEXT_TEMPLATE: str = '''
<div class="phylo-context">
    <h4>Phylogenetic Context</h4>
    <p class="phylo-intro">
        Novel taxa are placed within the taxonomic framework based on their nearest
        reference genomes. This provides context for understanding their evolutionary
        relationships.
    </p>
    <div class="phylo-list">
{context_items}
    </div>
</div>
'''

NOVEL_PHYLO_ITEM_TEMPLATE: str = '''<div class="phylo-item {type_class}">
    <span class="phylo-cluster">{cluster_id}</span>
    <span class="phylo-name">{suggested_name}</span>
    <span class="phylo-placement">{phylogenetic_context} (~{estimated_ani:.0f}% ANI to {closest_taxon})</span>
</div>'''

NOVEL_EMPTY_TEMPLATE: str = '''
<div class="novel-empty">
    <p class="empty-message">
        No novel diversity clusters found in this sample.
        This may indicate:
    </p>
    <ul>
        <li>All reads match known species in the reference database</li>
        <li>Novel reads did not meet the minimum cluster size threshold</li>
        <li>Sample has low novel diversity content</li>
    </ul>
</div>
'''

PHYLOGENETIC_HEATMAP_LEGEND_TEMPLATE: str = '''
<div class="phylo-heatmap-legend">
    <h4>Phylogenetic Context Heatmap Guide</h4>
    <div class="legend-content">
        <div class="legend-section">
            <span class="legend-title">Entity Types</span>
            <div class="legend-items">
                <div class="legend-item">
                    <span class="legend-marker ref-marker">GCF_*</span>
                    <span class="legend-desc">Reference genome from similarity matrix</span>
                </div>
                <div class="legend-item">
                    <span class="legend-marker novel-marker">[*] NSP_*</span>
                    <span class="legend-desc">Novel cluster (species or genus)</span>
                </div>
            </div>
        </div>
        <div class="legend-section">
            <span class="legend-title">ANI Thresholds (Nucleotide)</span>
            <div class="legend-items">
                <div class="legend-item">
                    <span class="legend-color" style="background: #d73027;"></span>
                    <span class="legend-desc">&gt;95% - Same species</span>
                </div>
                <div class="legend-item">
                    <span class="legend-color" style="background: #fdae61;"></span>
                    <span class="legend-desc">90-95% - Species boundary</span>
                </div>
                <div class="legend-item">
                    <span class="legend-color" style="background: #ffffbf;"></span>
                    <span class="legend-desc">80-90% - Same genus</span>
                </div>
                <div class="legend-item">
                    <span class="legend-color" style="background: #4575b4;"></span>
                    <span class="legend-desc">&lt;80% - Different genera</span>
                </div>
            </div>
        </div>
        <div class="legend-section">
            <span class="legend-title">AAI Thresholds (Protein)</span>
            <div class="legend-items">
                <div class="legend-item">
                    <span class="legend-color" style="background: #a50026;"></span>
                    <span class="legend-desc">&gt;80% - High similarity (same genus)</span>
                </div>
                <div class="legend-item">
                    <span class="legend-color" style="background: #fdae61;"></span>
                    <span class="legend-desc">65-80% - Genus boundary zone</span>
                </div>
                <div class="legend-item">
                    <span class="legend-color" style="background: #abd9e9;"></span>
                    <span class="legend-desc">50-65% - Same family</span>
                </div>
                <div class="legend-item">
                    <span class="legend-color" style="background: #313695;"></span>
                    <span class="legend-desc">&lt;50% - Different families</span>
                </div>
            </div>
        </div>
        <div class="legend-note">
            <strong>Note:</strong> Similarity values between novel clusters and references
            are estimated based on classification data. Values marked "(estimated)"
            in hover text are derived from the cluster's divergence from its
            nearest reference genome. AAI heatmaps are particularly useful for
            Novel Genus candidates where protein-level comparisons remain informative.
        </div>
    </div>
</div>
'''

PHYLOGENETIC_HEATMAP_INTRO_TEMPLATE: str = '''
<div class="phylo-heatmap-intro">
    <h4>Phylogenetic Placement of Novel Diversity</h4>
    <p>
        This heatmap shows novel clusters positioned alongside reference genomes
        based on estimated similarity relationships. Hierarchical clustering groups
        similar entities together, allowing you to see where novel taxa fall
        within the phylogenetic landscape of known species.
    </p>
    <div class="phylo-stats">
        <span class="stat-item"><strong>{n_references}</strong> reference genomes</span>
        <span class="stat-item"><strong>{n_clusters}</strong> novel clusters{clusters_note}</span>
    </div>
</div>
'''


# =============================================================================
# Builder Functions
# =============================================================================

def build_novel_summary_html(summary: NovelDiversitySummary) -> str:
    """
    Build the summary metrics HTML section.

    Args:
        summary: NovelDiversitySummary with aggregate statistics

    Returns:
        Formatted HTML string
    """
    return NOVEL_SUMMARY_TEMPLATE.format(
        total_clusters=summary.total_clusters,
        total_reads=summary.total_novel_reads,
        novel_species_clusters=summary.novel_species_clusters,
        novel_species_reads=summary.novel_species_reads,
        novel_genus_clusters=summary.novel_genus_clusters,
        novel_genus_reads=summary.novel_genus_reads,
        high_confidence=summary.high_confidence_clusters,
    )


def build_cluster_table_html(clusters: list[NovelCluster]) -> str:
    """
    Build the cluster details table HTML.

    Args:
        clusters: List of NovelCluster objects

    Returns:
        Formatted HTML string with table
    """
    if not clusters:
        return NOVEL_EMPTY_TEMPLATE

    rows = []
    for cluster in clusters:
        confidence_class = f"conf-{cluster.confidence.lower()}"
        type_class = "type-species" if cluster.taxonomic_call == "Novel Species" else "type-genus"

        # Build full taxonomy string for tooltip
        full_taxonomy = f"{cluster.nearest_species} ({cluster.nearest_genus}, {cluster.nearest_family})"

        row = NOVEL_CLUSTER_ROW_TEMPLATE.format(
            cluster_id=cluster.cluster_id,
            taxonomic_call=cluster.taxonomic_call,
            type_class=type_class,
            read_count=cluster.read_count,
            closest_taxon=cluster.closest_known_taxon,
            full_taxonomy=full_taxonomy,
            phylogenetic_context=cluster.phylogenetic_context,
            estimated_ani=cluster.estimated_ani,
            mean_uncertainty=cluster.mean_placement_uncertainty,
            confidence=cluster.confidence,
            confidence_class=confidence_class,
            suggested_name=cluster.suggested_name,
        )
        rows.append(row)

    return NOVEL_CLUSTER_TABLE_TEMPLATE.format(cluster_rows="\n".join(rows))


def build_phylogenetic_context_html(clusters: list[NovelCluster]) -> str:
    """
    Build the phylogenetic context section HTML.

    Args:
        clusters: List of NovelCluster objects

    Returns:
        Formatted HTML string
    """
    if not clusters:
        return ""

    items = []
    for cluster in clusters:
        type_class = "type-species" if cluster.taxonomic_call == "Novel Species" else "type-genus"

        item = NOVEL_PHYLO_ITEM_TEMPLATE.format(
            cluster_id=cluster.cluster_id,
            suggested_name=cluster.suggested_name,
            phylogenetic_context=cluster.phylogenetic_context,
            estimated_ani=cluster.estimated_ani,
            closest_taxon=cluster.closest_known_taxon,
            type_class=type_class,
        )
        items.append(item)

    return NOVEL_PHYLOGENETIC_CONTEXT_TEMPLATE.format(
        context_items="\n".join(items)
    )


def build_cluster_scatter_figure(
    clusters: list[NovelCluster],
    width: int = 800,
    height: int = 500,
) -> go.Figure:
    """
    Build a scatter plot of cluster quality.

    X-axis: Mean novelty index
    Y-axis: Mean discovery score (or confidence level if no discovery score)
    Marker size: Read count
    Color: Confidence rating

    Args:
        clusters: List of NovelCluster objects
        width: Plot width in pixels
        height: Plot height in pixels

    Returns:
        Plotly Figure object
    """
    import plotly.graph_objects as go

    if not clusters:
        fig = go.Figure()
        fig.update_layout(
            title="No clusters to display",
            template="plotly_white",
            width=width,
            height=height,
        )
        return fig

    # Group by confidence for coloring
    confidence_colors = {
        "High": "#22c55e",
        "Medium": "#f59e0b",
        "Low": "#94a3b8",
    }

    fig = go.Figure()

    for conf, color in confidence_colors.items():
        conf_clusters = [c for c in clusters if c.confidence == conf]
        if not conf_clusters:
            continue

        # Y-axis: use discovery score if available, else encode confidence as number
        if conf_clusters[0].mean_discovery_score is not None:
            y_vals = [c.mean_discovery_score for c in conf_clusters]
            y_title = "Mean Discovery Score"
        else:
            # Encode confidence as numeric
            conf_to_num = {"High": 80, "Medium": 50, "Low": 20}
            y_vals = [conf_to_num[c.confidence] for c in conf_clusters]
            y_title = "Confidence Level"

        # Scale marker sizes (min 10, max 50)
        max_reads = max(c.read_count for c in clusters)
        sizes = [
            10 + (c.read_count / max_reads) * 40
            for c in conf_clusters
        ]

        fig.add_trace(go.Scatter(
            x=[c.mean_novelty_index for c in conf_clusters],
            y=y_vals,
            mode="markers",
            name=conf,
            marker=dict(
                color=color,
                size=sizes,
                opacity=0.7,
                line=dict(width=1, color="white"),
            ),
            text=[
                f"{c.cluster_id}<br>{c.suggested_name}<br>{c.read_count} reads"
                for c in conf_clusters
            ],
            hovertemplate=(
                "<b>%{text}</b><br>"
                "Novelty: %{x:.1f}%<br>"
                f"{y_title}: %{{y:.1f}}<extra></extra>"
            ),
        ))

    fig.update_layout(
        title="Cluster Quality Overview",
        xaxis_title="Mean Novelty Index (%)",
        yaxis_title=y_title,
        template="plotly_white",
        width=width,
        height=height,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.2,
            xanchor="center",
            x=0.5,
        ),
        showlegend=True,
    )

    # Add novelty threshold lines
    fig.add_vline(
        x=5, line_dash="dot", line_color="#94a3b8",
        annotation_text="Novel Species (5%)",
        annotation_position="top",
    )
    fig.add_vline(
        x=20, line_dash="dot", line_color="#94a3b8",
        annotation_text="Novel Genus (20%)",
        annotation_position="top",
    )

    return fig


def build_sunburst_figure(
    clusters: list[NovelCluster],
    width: int = 600,
    height: int = 600,
) -> go.Figure:
    """
    Build a sunburst chart showing phylogenetic context.

    Center: Family level
    Middle ring: Genus level (known + "Novel" segments)
    Outer ring: Species/Novel clusters

    Args:
        clusters: List of NovelCluster objects
        width: Plot width in pixels
        height: Plot height in pixels

    Returns:
        Plotly Figure object
    """
    import plotly.graph_objects as go

    if not clusters:
        fig = go.Figure()
        fig.update_layout(
            title="No clusters to display",
            template="plotly_white",
            width=width,
            height=height,
        )
        return fig

    # Confidence to color mapping
    conf_colors = {
        "High": "#22c55e",
        "Medium": "#f59e0b",
        "Low": "#94a3b8",
    }

    # First pass: calculate totals for each family and genus
    family_totals: dict[str, int] = {}
    genus_totals: dict[str, int] = {}

    for cluster in clusters:
        family = cluster.nearest_family
        genus = cluster.nearest_genus
        genus_key = f"{family}/{genus}"

        family_totals[family] = family_totals.get(family, 0) + cluster.read_count
        genus_totals[genus_key] = genus_totals.get(genus_key, 0) + cluster.read_count

    # Build hierarchy: family -> genus -> cluster
    ids = []
    labels = []
    parents = []
    values = []
    colors = []
    custom_data = []

    # Track which families and genera have been added
    added_families: set[str] = set()
    added_genera: set[str] = set()

    # Add nodes in order: families first, then genera, then clusters
    for cluster in clusters:
        family = cluster.nearest_family
        genus = cluster.nearest_genus
        genus_key = f"{family}/{genus}"

        # Add family if not already added
        if family not in added_families:
            ids.append(family)
            labels.append(family)
            parents.append("")
            values.append(family_totals[family])
            colors.append("#e2e8f0")
            custom_data.append(f"Family: {family}")
            added_families.add(family)

        # Add genus if not already added
        if genus_key not in added_genera:
            ids.append(genus_key)
            labels.append(genus)
            parents.append(family)
            values.append(genus_totals[genus_key])
            colors.append("#cbd5e1")
            custom_data.append(f"Genus: {genus}")
            added_genera.add(genus_key)

    # Add cluster nodes (leaves)
    for cluster in clusters:
        family = cluster.nearest_family
        genus = cluster.nearest_genus
        genus_key = f"{family}/{genus}"
        cluster_node_id = f"{genus_key}/{cluster.cluster_id}"

        ids.append(cluster_node_id)
        labels.append(cluster.cluster_id)
        parents.append(genus_key)
        values.append(cluster.read_count)
        colors.append(conf_colors.get(cluster.confidence, "#94a3b8"))
        custom_data.append(f"{cluster.suggested_name} ({cluster.confidence})")

    fig = go.Figure(go.Sunburst(
        ids=ids,
        labels=labels,
        parents=parents,
        values=values,
        marker=dict(colors=colors),
        branchvalues="total",
        customdata=custom_data,
        hovertemplate=(
            "<b>%{label}</b><br>"
            "%{customdata}<br>"
            "Reads: %{value:,}<extra></extra>"
        ),
    ))

    fig.update_layout(
        title="Phylogenetic Context",
        template="plotly_white",
        width=width,
        height=height,
        margin=dict(t=50, l=10, r=10, b=10),
    )

    return fig


# =============================================================================
# CSS Styles for Novel Diversity Section
# =============================================================================

NOVEL_SECTION_CSS: str = '''
/* Novel Diversity Section Styles */
.novel-summary {
    margin-bottom: 2rem;
}

.novel-intro {
    margin-bottom: 1.5rem;
}

.novel-intro h3 {
    margin-bottom: 0.5rem;
    color: var(--text-primary);
}

.novel-intro p {
    color: var(--text-secondary);
    line-height: 1.6;
}

.novel-metrics {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
}

.novel-metric {
    background: var(--card-bg);
    border: 1px solid var(--border-color);
    border-radius: 8px;
    padding: 1.25rem;
    text-align: center;
}

.novel-metric.highlight {
    border-color: var(--accent-primary);
    background: linear-gradient(135deg, var(--card-bg), rgba(102, 126, 234, 0.1));
}

.novel-metric-value {
    font-size: 2rem;
    font-weight: 700;
    color: var(--text-primary);
}

.novel-metric.species .novel-metric-value { color: #f59e0b; }
.novel-metric.genus .novel-metric-value { color: #ef4444; }
.novel-metric.confidence .novel-metric-value { color: #22c55e; }

.novel-metric-label {
    font-size: 0.875rem;
    font-weight: 600;
    color: var(--text-secondary);
    margin-top: 0.25rem;
}

.novel-metric-hint {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.25rem;
}

/* Confidence Guide */
.novel-confidence-guide {
    margin: 1.5rem 0;
    padding: 1rem;
    background: var(--card-bg);
    border-radius: 8px;
    border: 1px solid var(--border-color);
}

.novel-confidence-guide h4 {
    margin-bottom: 1rem;
    color: var(--text-primary);
}

.confidence-grid {
    display: grid;
    gap: 0.75rem;
}

.confidence-item {
    display: grid;
    grid-template-columns: 80px 1fr 1fr;
    gap: 1rem;
    align-items: center;
    padding: 0.5rem;
    border-radius: 4px;
}

.confidence-item.high { background: rgba(34, 197, 94, 0.1); }
.confidence-item.medium { background: rgba(245, 158, 11, 0.1); }
.confidence-item.low { background: rgba(148, 163, 184, 0.1); }

.conf-badge {
    display: inline-block;
    padding: 0.25rem 0.75rem;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.75rem;
    text-transform: uppercase;
}

.confidence-item.high .conf-badge { background: #22c55e; color: white; }
.confidence-item.medium .conf-badge { background: #f59e0b; color: white; }
.confidence-item.low .conf-badge { background: #94a3b8; color: white; }

.conf-criteria {
    font-size: 0.875rem;
    color: var(--text-secondary);
}

.conf-action {
    font-size: 0.875rem;
    color: var(--text-muted);
    font-style: italic;
}

/* Cluster Table */
.novel-cluster-table {
    margin: 1.5rem 0;
}

.novel-cluster-table h4 {
    margin-bottom: 1rem;
}

.cluster-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.875rem;
}

.cluster-table th {
    text-align: left;
    padding: 0.75rem;
    background: var(--card-bg);
    border-bottom: 2px solid var(--border-color);
    font-weight: 600;
    color: var(--text-secondary);
}

.cluster-table td {
    padding: 0.75rem;
    border-bottom: 1px solid var(--border-color);
}

.cluster-row:hover {
    background: rgba(102, 126, 234, 0.05);
}

.cluster-id {
    font-family: monospace;
    font-weight: 600;
}

.cluster-type.type-species { color: #f59e0b; }
.cluster-type.type-genus { color: #ef4444; }

.cluster-placement {
    max-width: 300px;
}

.placement-name {
    font-weight: 600;
    font-style: italic;
    color: var(--text-primary);
    margin-bottom: 0.25rem;
}

.placement-context {
    font-size: 0.85rem;
    color: var(--text-secondary);
    margin-bottom: 0.25rem;
}

.placement-ref {
    font-size: 0.8rem;
    color: var(--text-muted);
}

.cluster-ani {
    font-weight: 600;
    color: var(--accent-primary);
}

.conf-tag {
    display: inline-block;
    padding: 0.2rem 0.5rem;
    border-radius: 3px;
    font-size: 0.75rem;
    font-weight: 600;
}

.conf-tag.conf-high { background: #22c55e; color: white; }
.conf-tag.conf-medium { background: #f59e0b; color: white; }
.conf-tag.conf-low { background: #94a3b8; color: white; }

.cluster-name {
    font-style: italic;
    color: var(--text-secondary);
}

/* Phylogenetic Context */
.phylo-context {
    margin: 1.5rem 0;
}

.phylo-intro {
    color: var(--text-secondary);
    margin-bottom: 1rem;
}

.phylo-list {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
}

.phylo-item {
    display: grid;
    grid-template-columns: 80px 200px 1fr;
    gap: 1rem;
    padding: 0.75rem;
    background: var(--card-bg);
    border-radius: 4px;
    border-left: 3px solid;
}

.phylo-item.type-species { border-left-color: #f59e0b; }
.phylo-item.type-genus { border-left-color: #ef4444; }

.phylo-cluster {
    font-family: monospace;
    font-weight: 600;
    color: var(--text-primary);
}

.phylo-name {
    font-style: italic;
    color: var(--text-secondary);
}

.phylo-placement {
    color: var(--text-muted);
}

/* Empty state */
.novel-empty {
    text-align: center;
    padding: 3rem;
    color: var(--text-muted);
}

.novel-empty ul {
    text-align: left;
    display: inline-block;
    margin-top: 1rem;
}

/* Phylogenetic Heatmap Styles */
.phylo-heatmap-intro {
    margin: 1.5rem 0;
    padding: 1rem;
    background: var(--card-bg);
    border-radius: 8px;
    border: 1px solid var(--border-color);
}

.phylo-heatmap-intro h4 {
    margin-bottom: 0.75rem;
    color: var(--text-primary);
}

.phylo-heatmap-intro p {
    color: var(--text-secondary);
    line-height: 1.6;
    margin-bottom: 1rem;
}

.phylo-stats {
    display: flex;
    gap: 2rem;
}

.stat-item {
    color: var(--text-secondary);
}

.stat-item strong {
    color: var(--accent-primary);
}

.phylo-heatmap-legend {
    margin: 1rem 0;
    padding: 1rem;
    background: var(--card-bg);
    border-radius: 8px;
    border: 1px solid var(--border-color);
}

.phylo-heatmap-legend h4 {
    margin-bottom: 1rem;
    font-size: 0.9rem;
    color: var(--text-primary);
}

.legend-content {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1.5rem;
}

@media (max-width: 768px) {
    .legend-content {
        grid-template-columns: 1fr;
    }
}

.legend-section {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
}

.legend-title {
    font-weight: 600;
    font-size: 0.85rem;
    color: var(--text-secondary);
    margin-bottom: 0.25rem;
}

.legend-items {
    display: flex;
    flex-direction: column;
    gap: 0.4rem;
}

.legend-item {
    display: flex;
    align-items: center;
    gap: 0.75rem;
}

.legend-marker {
    font-family: monospace;
    font-size: 0.8rem;
    padding: 0.2rem 0.4rem;
    border-radius: 3px;
    white-space: nowrap;
}

.legend-marker.ref-marker {
    background: #e2e8f0;
    color: #475569;
}

.legend-marker.novel-marker {
    background: #fef3c7;
    color: #b45309;
    font-weight: 600;
}

.legend-color {
    width: 20px;
    height: 14px;
    border-radius: 2px;
    flex-shrink: 0;
}

.legend-desc {
    font-size: 0.85rem;
    color: var(--text-secondary);
}

.legend-note {
    grid-column: 1 / -1;
    margin-top: 0.5rem;
    padding: 0.75rem;
    background: rgba(102, 126, 234, 0.05);
    border-radius: 4px;
    font-size: 0.85rem;
    color: var(--text-secondary);
    line-height: 1.5;
}

.legend-note strong {
    color: var(--text-primary);
}
'''
