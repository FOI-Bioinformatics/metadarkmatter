"""
CSS styles for HTML reports.

Provides both light and dark theme options with scientific,
publication-ready styling for the unified report.
"""

from __future__ import annotations

# =============================================================================
# Light Theme
# =============================================================================

LIGHT_THEME: str = """
:root {
    --bg-page:       #f0f2f5;
    --bg-card:       #ffffff;
    --bg-primary:    #ffffff;
    --bg-secondary:  #f0f2f5;
    --bg-tertiary:   #e4e7eb;
    --bg-collapsible: #f8f9fa;
    --text-primary:  #1a1a2e;
    --text-secondary: #6b7280;
    --text-muted:    #9ca3af;
    --border-color:  #e4e7eb;
    --border-light:  #e4e7eb;
    --header-bg:     #1a1a2e;
    --accent-color:  #667eea;
    --accent-hover:  #5a6fd6;
    --success-color: #22c55e;
    --warning-color: #f59e0b;
    --danger-color:  #ef4444;
    --info-color:    #94a3b8;
    --offtarget-color: #3498db;
    --shadow-sm: 0 1px 2px rgba(0,0,0,0.05);
    --shadow-md: 0 2px 4px rgba(0,0,0,0.06);
    --shadow-lg: 0 4px 8px rgba(0,0,0,0.08);
    --radius-sm: 4px;
    --radius-md: 6px;
    --radius-lg: 8px;
}

* {
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    background-color: var(--bg-page);
    color: var(--text-primary);
    line-height: 1.5;
    font-size: 13px;
}

/* Header */
.report-header {
    background: var(--header-bg);
    color: white;
    padding: 1.25rem 2rem 1rem;
    text-align: center;
}

.report-header h1 {
    font-size: 1.25rem;
    font-weight: 600;
    margin-bottom: 0.5rem;
}

.report-header .metadata {
    display: flex;
    justify-content: center;
    gap: 2rem;
    font-size: 0.875rem;
    opacity: 0.9;
}

.report-header .metadata span {
    display: flex;
    align-items: center;
    gap: 0.5rem;
}

/* Navigation */
.tab-navigation {
    display: flex;
    justify-content: center;
    flex-wrap: wrap;
    gap: 0;
    background: var(--bg-card);
    border-bottom: 1px solid var(--border-color);
    padding: 0 1.25rem;
    position: sticky;
    top: 0;
    z-index: 100;
}

.tab-btn {
    padding: 0.75rem 1.25rem;
    border: none;
    background: transparent;
    cursor: pointer;
    font-size: 0.8rem;
    font-weight: 500;
    color: var(--text-secondary);
    border-bottom: 2px solid transparent;
    transition: all 0.2s ease;
    border-radius: 0;
}

.tab-btn:hover {
    background: transparent;
    color: var(--text-primary);
    border-bottom-color: var(--border-color);
}

.tab-btn.active {
    background: transparent;
    color: var(--text-primary);
    font-weight: 700;
    border-bottom-color: var(--accent-color);
}

/* Main Content */
.report-content {
    max-width: 1400px;
    margin: 0 auto;
    padding: 1.25rem;
}

.tab-section {
    display: none;
    animation: fadeIn 0.3s ease;
}

.tab-section.active {
    display: block;
}

@keyframes fadeIn {
    from { opacity: 0; transform: translateY(10px); }
    to { opacity: 1; transform: translateY(0); }
}

.section-title {
    font-size: 0.85rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 1rem;
    color: var(--text-primary);
    padding-bottom: 0.5rem;
    border-bottom: 1px solid var(--border-color);
}

/* Metric Cards */
.metric-cards {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
    margin-bottom: 2rem;
}

.metric-card {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1.25rem;
    box-shadow: var(--shadow-sm);
    text-align: center;
    border: 1px solid var(--border-light);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
}

.metric-card:hover {
    transform: translateY(-2px);
    box-shadow: var(--shadow-md);
}

.metric-value {
    font-size: 2rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1.2;
}

.metric-value.success { color: var(--success-color); }
.metric-value.warning { color: var(--warning-color); }
.metric-value.danger { color: var(--danger-color); }
.metric-value.info { color: var(--info-color); }
.metric-value.muted { color: var(--text-muted); }

.metric-label {
    font-size: 0.75rem;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-top: 0.25rem;
}

.metric-subtext {
    font-size: 0.8rem;
    color: var(--text-muted);
    margin-top: 0.25rem;
}

/* KPI Strip */
.kpi-strip {
    display: flex;
    gap: 0.75rem;
    margin-bottom: 1rem;
    flex-wrap: wrap;
}

.kpi-card {
    flex: 1;
    min-width: 120px;
    background: var(--bg-card);
    border-radius: var(--radius-md);
    padding: 0.75rem 1rem;
    border-left: 3px solid var(--accent-color);
    box-shadow: var(--shadow-sm);
}

.kpi-card.known { border-left-color: var(--success-color); }
.kpi-card.novel { border-left-color: var(--danger-color); }
.kpi-card.uncertain { border-left-color: var(--info-color); }
.kpi-card.off-target { border-left-color: var(--offtarget-color); }
.kpi-card.accent { border-left-color: var(--accent-color); }

.kpi-value {
    font-size: 1.75rem;
    font-weight: 700;
    font-variant-numeric: tabular-nums;
    line-height: 1.2;
    color: var(--text-primary);
}

.kpi-label {
    font-size: 0.7rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: var(--text-secondary);
    margin-top: 0.15rem;
}

/* Collapsible Panels */
.collapsible-panel {
    margin-bottom: 0.75rem;
    border-radius: var(--radius-md);
    border: 1px solid var(--border-color);
    overflow: hidden;
}

.collapsible-header {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.75rem 1rem;
    background: var(--bg-collapsible);
    cursor: pointer;
    user-select: none;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--text-primary);
}

.collapsible-header:hover {
    background: var(--bg-tertiary);
}

.collapsible-chevron {
    transition: transform 0.2s ease;
    font-size: 0.75rem;
    color: var(--text-muted);
}

.collapsible-panel.open .collapsible-chevron {
    transform: rotate(90deg);
}

.collapsible-body {
    max-height: 0;
    overflow: hidden;
    transition: max-height 0.3s ease;
    background: var(--bg-card);
}

.collapsible-panel.open .collapsible-body {
    max-height: 5000px;
}

.collapsible-content {
    padding: 1rem;
}

/* Two-Column Row */
.two-col-row {
    display: flex;
    gap: 0.75rem;
    margin-bottom: 0.75rem;
}

.two-col-row .col-left,
.two-col-row .col-right {
    min-width: 0;
}

@media (max-width: 768px) {
    .two-col-row {
        flex-direction: column;
    }
}

/* Mini Table (for top species) */
.top-species-card {
    background: var(--bg-card);
    border-radius: var(--radius-md);
    padding: 1rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-color);
    height: 100%;
}

.card-section-title {
    font-size: 0.8rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: var(--text-secondary);
    margin-bottom: 0.75rem;
}

.mini-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.8rem;
}

.mini-table th {
    text-align: left;
    font-weight: 600;
    font-size: 0.7rem;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: var(--text-muted);
    padding: 0.25rem 0.5rem;
    border-bottom: 1px solid var(--border-color);
}

.mini-table td {
    padding: 0.3rem 0.5rem;
    border-bottom: 1px solid var(--border-light);
    color: var(--text-primary);
}

.mini-table .species-name {
    font-style: italic;
    max-width: 200px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}

.mini-table .num {
    text-align: right;
    font-variant-numeric: tabular-nums;
    white-space: nowrap;
}

/* Category Grid (compact 2-row) */
.category-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 0.5rem;
    margin-bottom: 0.75rem;
}

.cat-grid-cell {
    background: var(--bg-card);
    border-radius: var(--radius-sm);
    padding: 0.6rem 0.75rem;
    border-left: 3px solid var(--border-color);
    box-shadow: var(--shadow-sm);
}

.cat-grid-cell.known { border-left-color: var(--success-color); }
.cat-grid-cell.novel-species { border-left-color: var(--warning-color); }
.cat-grid-cell.novel-genus { border-left-color: var(--danger-color); }
.cat-grid-cell.species-boundary { border-left-color: var(--info-color); }
.cat-grid-cell.ambiguous { border-left-color: var(--info-color); }
.cat-grid-cell.conserved { border-left-color: var(--info-color); }
.cat-grid-cell.unclassified { border-left-color: var(--text-muted); }
.cat-grid-cell.off-target { border-left-color: var(--offtarget-color); }

.cat-grid-name {
    font-size: 0.7rem;
    font-weight: 600;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.3px;
}

.cat-grid-count {
    font-size: 1.1rem;
    font-weight: 700;
    font-variant-numeric: tabular-nums;
    color: var(--text-primary);
}

.cat-grid-pct {
    font-size: 0.7rem;
    color: var(--text-muted);
    font-variant-numeric: tabular-nums;
}

@media (max-width: 768px) {
    .category-grid {
        grid-template-columns: repeat(2, 1fr);
    }
}

/* Diversity Summary (Hero Section) */
.diversity-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 2rem;
    margin-bottom: 2rem;
    box-shadow: var(--shadow-md);
    border: 1px solid var(--border-light);
}

.diversity-cards {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1.5rem;
    margin-bottom: 2rem;
}

.diversity-card {
    padding: 1.5rem;
    border-radius: var(--radius-md);
    text-align: center;
    border: 2px solid transparent;
}

.diversity-card.known {
    background: rgba(46, 204, 113, 0.1);
    border-color: var(--success-color);
}

.diversity-card.novel {
    background: rgba(231, 76, 60, 0.1);
    border-color: var(--danger-color);
}

.diversity-card.uncertain {
    background: rgba(149, 165, 166, 0.1);
    border-color: var(--info-color);
}

.diversity-card.off-target {
    background: rgba(52, 152, 219, 0.1);
    border-color: #3498db;
}

.diversity-pct {
    font-size: 2.5rem;
    font-weight: 700;
    line-height: 1.1;
}

.diversity-card.known .diversity-pct { color: var(--success-color); }
.diversity-card.novel .diversity-pct { color: var(--danger-color); }
.diversity-card.uncertain .diversity-pct { color: var(--info-color); }
.diversity-card.off-target .diversity-pct { color: #3498db; }

.diversity-label {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0.5rem 0 0.25rem;
}

.diversity-count {
    font-size: 0.9rem;
    color: var(--text-secondary);
}

.diversity-desc {
    font-size: 0.8rem;
    color: var(--text-muted);
    margin-top: 0.5rem;
}

.diversity-bar {
    display: flex;
    height: 32px;
    border-radius: var(--radius-md);
    overflow: hidden;
    margin-bottom: 0.75rem;
}

.bar-segment {
    transition: width 0.3s ease;
    min-width: 2px;
}

.bar-segment.known { background: var(--success-color); }
.bar-segment.novel { background: var(--danger-color); }
.bar-segment.uncertain { background: var(--info-color); }
.bar-segment.off-target { background: var(--offtarget-color); }

.diversity-bar-legend {
    display: flex;
    justify-content: center;
    gap: 1.5rem;
    font-size: 0.75rem;
    color: var(--text-secondary);
    margin-bottom: 1rem;
}

.legend-item {
    display: flex;
    align-items: center;
    gap: 0.5rem;
}

.legend-dot {
    width: 12px;
    height: 12px;
    border-radius: 50%;
}

.legend-dot.known { background: var(--success-color); }
.legend-dot.novel { background: var(--danger-color); }
.legend-dot.uncertain { background: var(--info-color); }
.legend-dot.off-target { background: var(--offtarget-color); }

/* Key Findings Cards */
.key-findings {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 2rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.findings-title {
    font-size: 1.1rem;
    font-weight: 600;
    margin-bottom: 1rem;
    color: var(--text-primary);
}

.findings-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
    gap: 1rem;
}

.finding-card {
    display: flex;
    gap: 0.75rem;
    padding: 1rem;
    border-radius: var(--radius-md);
    background: var(--bg-secondary);
    border-left: 3px solid var(--border-color);
}

.finding-card.species { border-left-color: var(--success-color); }
.finding-card.novel { border-left-color: var(--danger-color); }
.finding-card.confidence { border-left-color: var(--accent-color); }
.finding-card.action { border-left-color: var(--warning-color); }

.finding-icon {
    font-size: 1.5rem;
    line-height: 1;
    flex-shrink: 0;
}

.finding-headline {
    font-weight: 600;
    font-size: 0.95rem;
    margin-bottom: 0.25rem;
    color: var(--text-primary);
}

.finding-detail {
    font-size: 0.85rem;
    color: var(--text-secondary);
    line-height: 1.4;
}

.finding-link {
    font-size: 0.8rem;
    color: var(--accent-color);
    text-decoration: none;
    margin-top: 0.25rem;
    display: inline-block;
}

/* Category Breakdown */
.category-breakdown {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 2rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.breakdown-title {
    font-size: 1.1rem;
    font-weight: 600;
    margin-bottom: 0.5rem;
    color: var(--text-primary);
}

.breakdown-explanation {
    font-size: 0.85rem;
    color: var(--text-secondary);
    margin-bottom: 1.5rem;
}

.breakdown-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1.5rem;
}

.breakdown-group {
    padding: 1rem;
    border-radius: var(--radius-md);
    border-left: 4px solid transparent;
}

.breakdown-group.known {
    background: rgba(46, 204, 113, 0.05);
    border-left-color: var(--success-color);
}

.breakdown-group.novel {
    background: rgba(231, 76, 60, 0.05);
    border-left-color: var(--danger-color);
}

.breakdown-group.uncertain {
    background: rgba(149, 165, 166, 0.05);
    border-left-color: var(--info-color);
}

.group-header {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 0.75rem;
}

.group-dot {
    width: 10px;
    height: 10px;
    border-radius: 50%;
}

.group-dot.known { background: var(--success-color); }
.group-dot.novel { background: var(--danger-color); }
.group-dot.uncertain { background: var(--info-color); }

.group-name {
    font-weight: 600;
    font-size: 0.9rem;
    color: var(--text-primary);
}

.category-item {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 0.4rem 0;
    border-bottom: 1px solid var(--border-light);
}

.category-item:last-of-type {
    border-bottom: none;
}

.cat-name {
    font-size: 0.85rem;
    color: var(--text-primary);
}

.cat-count {
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--text-primary);
}

.cat-pct {
    font-size: 0.8rem;
    color: var(--text-muted);
    min-width: 45px;
    text-align: right;
}

.category-desc {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.75rem;
    font-style: italic;
}

/* Distributions Summary */
.distributions-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.dist-intro {
    margin-bottom: 1.5rem;
}

.dist-intro h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.dist-intro p {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
}

.dist-metrics {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 1rem;
}

.dist-metric {
    text-align: center;
    padding: 1rem;
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
}

.dist-metric-value {
    font-size: 1.75rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1.2;
}

.dist-metric-label {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-top: 0.25rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.dist-metric-hint {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.25rem;
}

/* Scatter Plot Interpretation Guide */
.scatter-interpretation {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem 1.25rem;
    margin-bottom: 1rem;
    border: 1px solid var(--border-light);
}

.scatter-interpretation h4 {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.interpretation-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1rem;
}

.interp-region {
    display: flex;
    align-items: flex-start;
    gap: 0.75rem;
    padding: 0.75rem;
    border-radius: var(--radius-sm);
    background: var(--bg-secondary);
}

.region-marker {
    width: 16px;
    height: 16px;
    border-radius: 4px;
    flex-shrink: 0;
    margin-top: 2px;
}

.interp-region.known .region-marker { background: var(--success-color); }
.interp-region.novel .region-marker { background: var(--danger-color); }
.interp-region.uncertain .region-marker { background: var(--info-color); }

.region-text {
    display: flex;
    flex-direction: column;
    gap: 0.25rem;
}

.region-text strong {
    font-size: 0.8rem;
    color: var(--text-primary);
}

.region-text span {
    font-size: 0.75rem;
    color: var(--text-secondary);
    line-height: 1.4;
}

/* Genomes Summary */
.genomes-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.genomes-intro {
    margin-bottom: 1.5rem;
}

.genomes-intro h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.genomes-intro p {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
}

.genomes-metrics {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 1rem;
}

.genome-metric {
    text-align: center;
    padding: 1rem;
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
}

.genome-metric.highlight {
    background: linear-gradient(135deg, rgba(102, 126, 234, 0.1), rgba(118, 75, 162, 0.1));
    border: 1px solid var(--accent-color);
}

.genome-metric-value {
    font-size: 1.75rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1.2;
}

.genome-metric-label {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-top: 0.25rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.genome-metric-hint {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.25rem;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}

/* Genome Interpretation */
.genome-interpretation {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem 1.25rem;
    margin-bottom: 1.5rem;
    border: 1px solid var(--border-light);
}

.genome-interpretation h4 {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.interp-items {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
}

.interp-item {
    display: flex;
    gap: 0.5rem;
    font-size: 0.8rem;
}

.interp-label {
    font-weight: 600;
    color: var(--text-primary);
    min-width: 160px;
    flex-shrink: 0;
}

.interp-desc {
    color: var(--text-secondary);
}

/* Plot Containers */
.plot-row {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
    gap: 1.5rem;
    margin-bottom: 1.5rem;
}

.plot-container {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.plot-container.full-width {
    grid-column: 1 / -1;
}

.plot-title {
    font-size: 0.95rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid var(--border-light);
}

.plot-description {
    font-size: 0.8rem;
    color: var(--text-secondary);
    margin-bottom: 0.75rem;
}

/* Data Tables */
.table-container {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
    overflow: hidden;
}

.table-controls {
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    margin-bottom: 1rem;
    padding-bottom: 1rem;
    border-bottom: 1px solid var(--border-light);
}

.table-search {
    flex: 1;
    min-width: 200px;
    padding: 0.5rem 0.75rem;
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    transition: border-color 0.2s ease;
}

.table-search:focus {
    outline: none;
    border-color: var(--accent-color);
    box-shadow: 0 0 0 3px rgba(13, 110, 253, 0.15);
}

.table-filter {
    padding: 0.5rem 0.75rem;
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    background: var(--bg-primary);
    cursor: pointer;
}

.export-btn {
    padding: 0.5rem 1rem;
    background: var(--accent-color);
    color: white;
    border: none;
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    font-weight: 500;
    cursor: pointer;
    transition: background 0.2s ease;
}

.export-btn:hover {
    background: var(--accent-hover);
}

.column-selector {
    position: relative;
}

.column-btn {
    padding: 0.5rem 1rem;
    background: var(--bg-tertiary);
    color: var(--text-primary);
    border: 1px solid var(--border-light);
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    font-weight: 500;
    cursor: pointer;
    transition: background 0.2s ease;
}

.column-btn:hover {
    background: var(--bg-secondary);
}

.column-menu {
    display: none;
    position: absolute;
    top: 100%;
    right: 0;
    margin-top: 0.5rem;
    padding: 0.75rem;
    background: var(--bg-primary);
    border: 1px solid var(--border-light);
    border-radius: var(--radius-md);
    box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    z-index: 100;
    min-width: 180px;
}

.column-menu.show {
    display: block;
}

.column-menu label {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.4rem 0.25rem;
    font-size: 0.85rem;
    cursor: pointer;
    white-space: nowrap;
}

.column-menu label:hover {
    background: var(--bg-secondary);
    border-radius: var(--radius-sm);
}

.column-menu input[type="checkbox"] {
    width: 16px;
    height: 16px;
    cursor: pointer;
}

.data-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.8rem;
}

.data-table th,
.data-table td {
    padding: 0.6rem 0.75rem;
    text-align: left;
    border-bottom: 1px solid var(--border-light);
}

.data-table th {
    background: var(--bg-secondary);
    font-weight: 600;
    color: var(--text-primary);
    cursor: pointer;
    user-select: none;
    white-space: nowrap;
}

.data-table th:hover {
    background: var(--bg-tertiary);
}

.data-table th .sort-icon {
    margin-left: 0.25rem;
    opacity: 0.5;
}

.data-table th.sorted .sort-icon {
    opacity: 1;
}

.data-table tbody tr:hover {
    background: var(--bg-secondary);
}

.data-table .cell-known { color: var(--success-color); font-weight: 500; }
.data-table .cell-novel-species { color: var(--warning-color); font-weight: 500; }
.data-table .cell-novel-genus { color: var(--danger-color); font-weight: 500; }
.data-table .cell-conserved { color: var(--info-color); }

.table-pagination {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding-top: 1rem;
    border-top: 1px solid var(--border-light);
    margin-top: 1rem;
}

.page-info {
    font-size: 0.8rem;
    color: var(--text-secondary);
}

.page-buttons {
    display: flex;
    gap: 0.5rem;
}

.page-btn {
    padding: 0.4rem 0.75rem;
    border: 1px solid var(--border-color);
    background: var(--bg-primary);
    border-radius: var(--radius-sm);
    font-size: 0.8rem;
    cursor: pointer;
    transition: all 0.2s ease;
}

.page-btn:hover:not(:disabled) {
    background: var(--bg-secondary);
    border-color: var(--accent-color);
}

.page-btn:disabled {
    opacity: 0.5;
    cursor: not-allowed;
}

.page-number {
    font-size: 0.8rem;
    color: var(--text-secondary);
    padding: 0 0.5rem;
}

/* Data Tab Summary */
.data-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1.25rem;
    margin-bottom: 1.25rem;
    border: 1px solid var(--border-light);
}

.data-intro h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 0.5rem 0;
}

.data-intro p {
    font-size: 0.875rem;
    color: var(--text-secondary);
    margin: 0;
    line-height: 1.5;
}

.data-stats {
    display: flex;
    gap: 1.5rem;
    margin-top: 1rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border-light);
}

.data-stat {
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 0.5rem 1rem;
    border-radius: var(--radius-sm);
    background: var(--bg-secondary);
}

.data-stat .stat-value {
    font-size: 1.25rem;
    font-weight: 700;
    color: var(--text-primary);
}

.data-stat .stat-label {
    font-size: 0.75rem;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.data-stat.known { border-left: 3px solid var(--success-color); }
.data-stat.novel { border-left: 3px solid var(--warning-color); }
.data-stat.uncertain { border-left: 3px solid var(--text-muted); }
.data-stat.off-target { border-left: 3px solid #3498db; }

/* Quick Filters */
.quick-filters {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 1rem;
    padding: 0.75rem 1rem;
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
}

.filter-label {
    font-size: 0.8rem;
    color: var(--text-secondary);
    font-weight: 500;
    margin-right: 0.5rem;
}

.filter-chip {
    padding: 0.35rem 0.75rem;
    border: 1px solid var(--border-color);
    border-radius: 20px;
    background: var(--bg-primary);
    font-size: 0.75rem;
    cursor: pointer;
    transition: all 0.2s ease;
    white-space: nowrap;
}

.filter-chip:hover {
    border-color: var(--accent-color);
    background: var(--bg-tertiary);
}

.filter-chip.active {
    background: var(--accent-color);
    color: white;
    border-color: var(--accent-color);
}

.filter-chip.known.active { background: var(--success-color); border-color: var(--success-color); }
.filter-chip.novel-species.active { background: var(--warning-color); border-color: var(--warning-color); }
.filter-chip.novel-genus.active { background: var(--danger-color); border-color: var(--danger-color); }
.filter-chip.species-boundary.active { background: #8b5cf6; border-color: #8b5cf6; }
.filter-chip.ambiguous.active { background: var(--text-muted); border-color: var(--text-muted); }
.filter-chip.ambiguous-wg.active { background: #6b7280; border-color: #6b7280; }
.filter-chip.conserved.active { background: var(--info-color); border-color: var(--info-color); }
.filter-chip.unclassified.active { background: #d1d5db; border-color: #d1d5db; color: #1a1a2e; }
.filter-chip.off-target.active { background: #3498db; border-color: #3498db; }
.filter-chip.range.active { background: var(--accent-color); border-color: var(--accent-color); }

/* Filter Rows */
.quick-filters {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
    margin-bottom: 1rem;
    padding: 1rem;
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    border: 1px solid var(--border-light);
}

.filter-row {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    flex-wrap: wrap;
}

.filter-row .filter-label {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--text-secondary);
    min-width: 80px;
}

/* Truncation Notice */
.truncation-notice {
    background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%);
    border: 1px solid #f59e0b;
    border-radius: var(--radius-md);
    padding: 0.75rem 1rem;
    margin-bottom: 1rem;
    font-size: 0.85rem;
    color: #92400e;
}

.truncation-notice strong {
    color: #78350f;
}

/* Column Guide */
.column-guide {
    margin-bottom: 1rem;
    background: var(--bg-primary);
    border: 1px solid var(--border-light);
    border-radius: var(--radius-md);
}

.column-guide summary {
    padding: 0.75rem 1rem;
    font-size: 0.85rem;
    font-weight: 500;
    color: var(--text-secondary);
    cursor: pointer;
    user-select: none;
}

.column-guide summary:hover {
    background: var(--bg-secondary);
}

.column-guide[open] summary {
    border-bottom: 1px solid var(--border-light);
}

.guide-content {
    padding: 0.75rem 1rem;
}

.guide-item {
    display: flex;
    gap: 1rem;
    padding: 0.5rem 0;
    font-size: 0.8rem;
    border-bottom: 1px solid var(--border-light);
}

.guide-item:last-child {
    border-bottom: none;
}

.guide-col {
    font-weight: 600;
    color: var(--text-primary);
    min-width: 100px;
}

.guide-desc {
    color: var(--text-secondary);
}

.guide-section {
    margin-bottom: 1rem;
}

.guide-section:last-child {
    margin-bottom: 0;
}

.guide-section h4 {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--accent-color);
    text-transform: uppercase;
    letter-spacing: 0.03em;
    margin-bottom: 0.5rem;
    padding-bottom: 0.25rem;
    border-bottom: 1px solid var(--border-light);
}

.guide-section.enhanced-guide {
    background: var(--bg-secondary);
    padding: 0.75rem;
    border-radius: var(--radius-sm);
    margin-top: 0.5rem;
}

.guide-section.enhanced-guide h4 {
    color: var(--warning-color);
}

/* Enhanced Table Controls */
.controls-row {
    display: flex;
    gap: 0.75rem;
    width: 100%;
    align-items: center;
}

.controls-right {
    display: flex;
    gap: 0.5rem;
    margin-left: auto;
}

.filter-status {
    display: flex;
    align-items: center;
    gap: 1rem;
    font-size: 0.85rem;
}

.results-count {
    color: var(--text-secondary);
}

.results-count strong {
    color: var(--text-primary);
    font-weight: 600;
}

.clear-btn {
    padding: 0.3rem 0.6rem;
    background: transparent;
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    font-size: 0.75rem;
    color: var(--text-secondary);
    cursor: pointer;
    transition: all 0.2s ease;
}

.clear-btn:hover {
    background: var(--danger-color);
    color: white;
    border-color: var(--danger-color);
}

.menu-header {
    font-size: 0.75rem;
    font-weight: 600;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.03em;
    padding: 0.25rem 0.25rem 0.5rem 0.25rem;
    border-bottom: 1px solid var(--border-light);
    margin-bottom: 0.5rem;
}

.table-wrapper {
    overflow-x: auto;
}

/* Methods Footer Panel */
.methods-footer {
    max-width: 1400px;
    margin: 1rem auto 0;
    padding: 0 1.25rem 1.25rem;
}

/* Footer */
.report-footer {
    text-align: center;
    padding: 2rem;
    color: var(--text-muted);
    font-size: 0.8rem;
    border-top: 1px solid var(--border-color);
    margin-top: 2rem;
}

.report-footer a {
    color: var(--accent-color);
    text-decoration: none;
}

/* Responsive */
@media (max-width: 768px) {
    .report-header {
        padding: 1.5rem 1rem;
    }

    .report-header .metadata {
        flex-direction: column;
        gap: 0.5rem;
    }

    .tab-navigation {
        padding: 0.5rem;
    }

    .tab-btn {
        padding: 0.4rem 0.75rem;
        font-size: 0.8rem;
    }

    .report-content {
        padding: 1rem;
    }

    .plot-row {
        grid-template-columns: 1fr;
    }

    .metric-cards {
        grid-template-columns: repeat(2, 1fr);
    }

    .diversity-summary {
        padding: 1.25rem;
    }

    .diversity-cards {
        grid-template-columns: 1fr;
        gap: 1rem;
    }

    .diversity-pct {
        font-size: 2rem;
    }

    .breakdown-grid {
        grid-template-columns: 1fr;
        gap: 1rem;
    }

    .diversity-bar-legend {
        flex-wrap: wrap;
        gap: 1rem;
    }

    .dist-metrics {
        grid-template-columns: repeat(2, 1fr);
    }

    .interpretation-grid {
        grid-template-columns: 1fr;
    }

    .genomes-metrics {
        grid-template-columns: repeat(2, 1fr);
    }

    .interp-item {
        flex-direction: column;
        gap: 0.25rem;
    }

    .interp-label {
        min-width: auto;
    }

    .metric-grid {
        grid-template-columns: repeat(2, 1fr);
    }

    .confidence-cards {
        grid-template-columns: 1fr;
    }

    .type-cards {
        grid-template-columns: 1fr;
    }
}

/* Enhanced Scoring Section */
.enhanced-scoring-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.enhanced-scoring-summary .summary-intro h3 {
    font-size: 1.2rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.enhanced-scoring-summary .intro-text {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
    margin-bottom: 1.5rem;
}

.metric-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 1rem;
}

.metric-box {
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
    padding: 1rem;
    text-align: center;
}

.metric-box.highlight {
    background: linear-gradient(135deg, #667eea22 0%, #764ba222 100%);
    border: 1px solid #667eea44;
}

.metric-box .metric-title {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.03em;
    margin-bottom: 0.5rem;
}

.metric-box .metric-value {
    font-size: 2rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1.2;
}

.metric-box .metric-detail {
    font-size: 0.85rem;
    color: var(--text-primary);
    margin-top: 0.25rem;
}

.metric-box .metric-note {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.5rem;
    line-height: 1.4;
}

/* Confidence Dimensions */
.confidence-dimensions {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.confidence-dimensions h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.confidence-dimensions .section-intro {
    font-size: 0.9rem;
    color: var(--text-secondary);
    margin-bottom: 1.5rem;
    line-height: 1.5;
}

.confidence-cards {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1rem;
}

.confidence-card {
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
    padding: 1.25rem;
    text-align: center;
}

.confidence-card .card-header {
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.confidence-card .card-value {
    font-size: 2.5rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1;
    margin-bottom: 0.75rem;
}

.confidence-card .card-desc {
    font-size: 0.8rem;
    color: var(--text-secondary);
    line-height: 1.4;
    margin-bottom: 0.75rem;
}

.confidence-card .card-scale {
    display: flex;
    justify-content: space-between;
    font-size: 0.7rem;
    padding-top: 0.5rem;
    border-top: 1px solid var(--border-light);
}

.confidence-card .card-scale .low { color: var(--danger-color); }
.confidence-card .card-scale .mid { color: var(--warning-color); }
.confidence-card .card-scale .high { color: var(--success-color); }

/* Discovery Score Guide */
.discovery-guide {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.discovery-guide h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 1rem;
}

.guide-table {
    width: 100%;
    border-collapse: collapse;
}

.guide-table th,
.guide-table td {
    padding: 0.75rem 1rem;
    text-align: left;
    border-bottom: 1px solid var(--border-light);
}

.guide-table th {
    font-weight: 600;
    font-size: 0.8rem;
    text-transform: uppercase;
    color: var(--text-secondary);
    background: var(--bg-secondary);
}

.guide-table tr.priority-high td { border-left: 4px solid var(--success-color); }
.guide-table tr.priority-medium td { border-left: 4px solid var(--warning-color); }
.guide-table tr.priority-low td { border-left: 4px solid #f59e0b66; }
.guide-table tr.priority-uncertain td { border-left: 4px solid var(--text-muted); }

/* Uncertainty Types */
.uncertainty-types {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.uncertainty-types h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 1rem;
}

.type-cards {
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 1rem;
    margin-bottom: 1rem;
}

.type-card {
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
    padding: 1.25rem;
}

.type-card.measured {
    border-left: 4px solid var(--success-color);
}

.type-card.inferred {
    border-left: 4px solid var(--warning-color);
}

.type-card .type-header {
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.type-card .type-count {
    font-size: 1.5rem;
    font-weight: 700;
    color: var(--accent-color);
    margin-bottom: 0.5rem;
}

.type-card .type-desc {
    font-size: 0.85rem;
    color: var(--text-secondary);
    line-height: 1.4;
}

.uncertainty-note {
    font-size: 0.85rem;
    color: var(--text-secondary);
    background: var(--bg-secondary);
    padding: 1rem;
    border-radius: var(--radius-md);
    line-height: 1.5;
}

/* Print styles */
@media print {
    .tab-navigation {
        display: none;
    }

    .tab-section {
        display: block !important;
        page-break-inside: avoid;
    }

    .report-header {
        background: #1a1a2e;
        -webkit-print-color-adjust: exact;
        print-color-adjust: exact;
    }
}
"""

# =============================================================================
# Dark Theme
# =============================================================================

DARK_THEME: str = """
:root {
    --bg-page:       #0f0f1a;
    --bg-card:       #1e1e2e;
    --bg-primary:    #1a1a2e;
    --bg-secondary:  #16213e;
    --bg-tertiary:   #0f3460;
    --bg-collapsible: #252538;
    --text-primary:  #eaeaea;
    --text-secondary: #94a3b8;
    --text-muted:    #64748b;
    --border-color:  #334155;
    --border-light:  #1e293b;
    --header-bg:     #0f0f1a;
    --accent-color:  #818cf8;
    --accent-hover:  #6366f1;
    --success-color: #34d399;
    --warning-color: #fbbf24;
    --danger-color:  #f87171;
    --info-color:    #94a3b8;
    --offtarget-color: #3498db;
    --shadow-sm: 0 1px 2px rgba(0,0,0,0.3);
    --shadow-md: 0 2px 4px rgba(0,0,0,0.4);
    --shadow-lg: 0 4px 8px rgba(0,0,0,0.5);
    --radius-sm: 4px;
    --radius-md: 6px;
    --radius-lg: 8px;
}

* {
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    background-color: var(--bg-page);
    color: var(--text-primary);
    line-height: 1.5;
    font-size: 13px;
}

.report-header {
    background: var(--header-bg);
    color: white;
    padding: 1.25rem 2rem 1rem;
    text-align: center;
}

.report-header h1 {
    font-size: 1.25rem;
    font-weight: 600;
    margin-bottom: 0.5rem;
}

.report-header .metadata {
    display: flex;
    justify-content: center;
    gap: 2rem;
    font-size: 0.875rem;
    opacity: 0.9;
}

.tab-navigation {
    display: flex;
    justify-content: center;
    flex-wrap: wrap;
    gap: 0;
    background: var(--bg-card);
    border-bottom: 1px solid var(--border-color);
    padding: 0 1.25rem;
    position: sticky;
    top: 0;
    z-index: 100;
}

.tab-btn {
    padding: 0.75rem 1.25rem;
    border: none;
    background: transparent;
    cursor: pointer;
    font-size: 0.8rem;
    font-weight: 500;
    color: var(--text-secondary);
    border-bottom: 2px solid transparent;
    transition: all 0.2s ease;
    border-radius: 0;
}

.tab-btn:hover {
    background: transparent;
    color: var(--text-primary);
    border-bottom-color: var(--border-color);
}

.tab-btn.active {
    background: transparent;
    color: var(--text-primary);
    font-weight: 700;
    border-bottom-color: var(--accent-color);
}

.report-content {
    max-width: 1400px;
    margin: 0 auto;
    padding: 1.25rem;
}

.tab-section {
    display: none;
    animation: fadeIn 0.3s ease;
}

.tab-section.active {
    display: block;
}

@keyframes fadeIn {
    from { opacity: 0; transform: translateY(10px); }
    to { opacity: 1; transform: translateY(0); }
}

.section-title {
    font-size: 0.85rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 1rem;
    color: var(--text-primary);
    padding-bottom: 0.5rem;
    border-bottom: 1px solid var(--border-color);
}

.metric-cards {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
    margin-bottom: 2rem;
}

.metric-card {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1.25rem;
    box-shadow: var(--shadow-sm);
    text-align: center;
    border: 1px solid var(--border-color);
}

.metric-value {
    font-size: 2rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1.2;
}

.metric-value.success { color: var(--success-color); }
.metric-value.warning { color: var(--warning-color); }
.metric-value.danger { color: var(--danger-color); }
.metric-value.info { color: var(--info-color); }
.metric-value.muted { color: var(--text-muted); }

.metric-label {
    font-size: 0.75rem;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-top: 0.25rem;
}

/* KPI Strip - Dark Theme */
.kpi-strip {
    display: flex;
    gap: 0.75rem;
    margin-bottom: 1rem;
    flex-wrap: wrap;
}

.kpi-card {
    flex: 1;
    min-width: 120px;
    background: var(--bg-card);
    border-radius: var(--radius-md);
    padding: 0.75rem 1rem;
    border-left: 3px solid var(--accent-color);
    box-shadow: var(--shadow-sm);
}

.kpi-card.known { border-left-color: var(--success-color); }
.kpi-card.novel { border-left-color: var(--danger-color); }
.kpi-card.uncertain { border-left-color: var(--info-color); }
.kpi-card.off-target { border-left-color: var(--offtarget-color); }
.kpi-card.accent { border-left-color: var(--accent-color); }

.kpi-value {
    font-size: 1.75rem;
    font-weight: 700;
    font-variant-numeric: tabular-nums;
    line-height: 1.2;
    color: var(--text-primary);
}

.kpi-label {
    font-size: 0.7rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: var(--text-secondary);
    margin-top: 0.15rem;
}

/* Collapsible Panels - Dark Theme */
.collapsible-panel {
    margin-bottom: 0.75rem;
    border-radius: var(--radius-md);
    border: 1px solid var(--border-color);
    overflow: hidden;
}

.collapsible-header {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.75rem 1rem;
    background: var(--bg-collapsible);
    cursor: pointer;
    user-select: none;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--text-primary);
}

.collapsible-header:hover {
    background: var(--bg-tertiary);
}

.collapsible-chevron {
    transition: transform 0.2s ease;
    font-size: 0.75rem;
    color: var(--text-muted);
}

.collapsible-panel.open .collapsible-chevron {
    transform: rotate(90deg);
}

.collapsible-body {
    max-height: 0;
    overflow: hidden;
    transition: max-height 0.3s ease;
    background: var(--bg-card);
}

.collapsible-panel.open .collapsible-body {
    max-height: 5000px;
}

.collapsible-content {
    padding: 1rem;
}

/* Two-Column Row - Dark Theme */
.two-col-row {
    display: flex;
    gap: 0.75rem;
    margin-bottom: 0.75rem;
}

.two-col-row .col-left,
.two-col-row .col-right {
    min-width: 0;
}

@media (max-width: 768px) {
    .two-col-row {
        flex-direction: column;
    }
}

/* Mini Table - Dark Theme */
.top-species-card {
    background: var(--bg-card);
    border-radius: var(--radius-md);
    padding: 1rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-color);
    height: 100%;
}

.card-section-title {
    font-size: 0.8rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: var(--text-secondary);
    margin-bottom: 0.75rem;
}

.mini-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.8rem;
}

.mini-table th {
    text-align: left;
    font-weight: 600;
    font-size: 0.7rem;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: var(--text-muted);
    padding: 0.25rem 0.5rem;
    border-bottom: 1px solid var(--border-color);
}

.mini-table td {
    padding: 0.3rem 0.5rem;
    border-bottom: 1px solid var(--border-light);
    color: var(--text-primary);
}

.mini-table .species-name {
    font-style: italic;
    max-width: 200px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}

.mini-table .num {
    text-align: right;
    font-variant-numeric: tabular-nums;
    white-space: nowrap;
}

/* Category Grid - Dark Theme */
.category-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 0.5rem;
    margin-bottom: 0.75rem;
}

.cat-grid-cell {
    background: var(--bg-card);
    border-radius: var(--radius-sm);
    padding: 0.6rem 0.75rem;
    border-left: 3px solid var(--border-color);
    box-shadow: var(--shadow-sm);
}

.cat-grid-cell.known { border-left-color: var(--success-color); }
.cat-grid-cell.novel-species { border-left-color: var(--warning-color); }
.cat-grid-cell.novel-genus { border-left-color: var(--danger-color); }
.cat-grid-cell.species-boundary { border-left-color: var(--info-color); }
.cat-grid-cell.ambiguous { border-left-color: var(--info-color); }
.cat-grid-cell.conserved { border-left-color: var(--info-color); }
.cat-grid-cell.unclassified { border-left-color: var(--text-muted); }
.cat-grid-cell.off-target { border-left-color: var(--offtarget-color); }

.cat-grid-name {
    font-size: 0.7rem;
    font-weight: 600;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.3px;
}

.cat-grid-count {
    font-size: 1.1rem;
    font-weight: 700;
    font-variant-numeric: tabular-nums;
    color: var(--text-primary);
}

.cat-grid-pct {
    font-size: 0.7rem;
    color: var(--text-muted);
    font-variant-numeric: tabular-nums;
}

@media (max-width: 768px) {
    .category-grid {
        grid-template-columns: repeat(2, 1fr);
    }
}

/* Diversity Summary (Hero Section) - Dark Theme */
.diversity-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 2rem;
    margin-bottom: 2rem;
    box-shadow: var(--shadow-md);
    border: 1px solid var(--border-color);
}

.diversity-cards {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1.5rem;
    margin-bottom: 2rem;
}

.diversity-card {
    padding: 1.5rem;
    border-radius: var(--radius-md);
    text-align: center;
    border: 2px solid transparent;
}

.diversity-card.known {
    background: rgba(52, 211, 153, 0.1);
    border-color: var(--success-color);
}

.diversity-card.novel {
    background: rgba(248, 113, 113, 0.1);
    border-color: var(--danger-color);
}

.diversity-card.uncertain {
    background: rgba(148, 163, 184, 0.1);
    border-color: var(--info-color);
}

.diversity-card.off-target {
    background: rgba(52, 152, 219, 0.1);
    border-color: #3498db;
}

.diversity-pct {
    font-size: 2.5rem;
    font-weight: 700;
    line-height: 1.1;
}

.diversity-card.known .diversity-pct { color: var(--success-color); }
.diversity-card.novel .diversity-pct { color: var(--danger-color); }
.diversity-card.uncertain .diversity-pct { color: var(--info-color); }
.diversity-card.off-target .diversity-pct { color: #3498db; }

.diversity-label {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0.5rem 0 0.25rem;
}

.diversity-count {
    font-size: 0.9rem;
    color: var(--text-secondary);
}

.diversity-desc {
    font-size: 0.8rem;
    color: var(--text-muted);
    margin-top: 0.5rem;
}

.diversity-bar {
    display: flex;
    height: 32px;
    border-radius: var(--radius-md);
    overflow: hidden;
    margin-bottom: 0.75rem;
}

.bar-segment {
    transition: width 0.3s ease;
    min-width: 2px;
}

.bar-segment.known { background: var(--success-color); }
.bar-segment.novel { background: var(--danger-color); }
.bar-segment.uncertain { background: var(--info-color); }
.bar-segment.off-target { background: var(--offtarget-color); }

.diversity-bar-legend {
    display: flex;
    justify-content: center;
    gap: 1.5rem;
    font-size: 0.75rem;
    color: var(--text-secondary);
    margin-bottom: 1rem;
}

.legend-item {
    display: flex;
    align-items: center;
    gap: 0.5rem;
}

.legend-dot {
    width: 12px;
    height: 12px;
    border-radius: 50%;
}

.legend-dot.known { background: var(--success-color); }
.legend-dot.novel { background: var(--danger-color); }
.legend-dot.uncertain { background: var(--info-color); }
.legend-dot.off-target { background: var(--offtarget-color); }

/* Category Breakdown - Dark Theme */
.category-breakdown {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 2rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-color);
}

.breakdown-title {
    font-size: 1.1rem;
    font-weight: 600;
    margin-bottom: 0.5rem;
    color: var(--text-primary);
}

.breakdown-explanation {
    font-size: 0.85rem;
    color: var(--text-secondary);
    margin-bottom: 1.5rem;
}

.breakdown-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1.5rem;
}

.breakdown-group {
    padding: 1rem;
    border-radius: var(--radius-md);
    border-left: 4px solid transparent;
}

.breakdown-group.known {
    background: rgba(52, 211, 153, 0.05);
    border-left-color: var(--success-color);
}

.breakdown-group.novel {
    background: rgba(248, 113, 113, 0.05);
    border-left-color: var(--danger-color);
}

.breakdown-group.uncertain {
    background: rgba(148, 163, 184, 0.05);
    border-left-color: var(--info-color);
}

.group-header {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 0.75rem;
}

.group-dot {
    width: 10px;
    height: 10px;
    border-radius: 50%;
}

.group-dot.known { background: var(--success-color); }
.group-dot.novel { background: var(--danger-color); }
.group-dot.uncertain { background: var(--info-color); }

.group-name {
    font-weight: 600;
    font-size: 0.9rem;
    color: var(--text-primary);
}

.category-item {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 0.4rem 0;
    border-bottom: 1px solid var(--border-color);
}

.category-item:last-of-type {
    border-bottom: none;
}

.cat-name {
    font-size: 0.85rem;
    color: var(--text-primary);
}

.cat-count {
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--text-primary);
}

.cat-pct {
    font-size: 0.8rem;
    color: var(--text-muted);
    min-width: 45px;
    text-align: right;
}

.category-desc {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.75rem;
    font-style: italic;
}

/* Distributions Summary - Dark Theme */
.distributions-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-color);
}

.dist-intro {
    margin-bottom: 1.5rem;
}

.dist-intro h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.dist-intro p {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
}

.dist-metrics {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 1rem;
}

.dist-metric {
    text-align: center;
    padding: 1rem;
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
}

.dist-metric-value {
    font-size: 1.75rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1.2;
}

.dist-metric-label {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-top: 0.25rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.dist-metric-hint {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.25rem;
}

/* Scatter Plot Interpretation Guide - Dark Theme */
.scatter-interpretation {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem 1.25rem;
    margin-bottom: 1rem;
    border: 1px solid var(--border-color);
}

.scatter-interpretation h4 {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.interpretation-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1rem;
}

.interp-region {
    display: flex;
    align-items: flex-start;
    gap: 0.75rem;
    padding: 0.75rem;
    border-radius: var(--radius-sm);
    background: var(--bg-secondary);
}

.region-marker {
    width: 16px;
    height: 16px;
    border-radius: 4px;
    flex-shrink: 0;
    margin-top: 2px;
}

.interp-region.known .region-marker { background: var(--success-color); }
.interp-region.novel .region-marker { background: var(--danger-color); }
.interp-region.uncertain .region-marker { background: var(--info-color); }

.region-text {
    display: flex;
    flex-direction: column;
    gap: 0.25rem;
}

.region-text strong {
    font-size: 0.8rem;
    color: var(--text-primary);
}

.region-text span {
    font-size: 0.75rem;
    color: var(--text-secondary);
    line-height: 1.4;
}

/* Genomes Summary - Dark Theme */
.genomes-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-color);
}

.genomes-intro {
    margin-bottom: 1.5rem;
}

.genomes-intro h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.genomes-intro p {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
}

.genomes-metrics {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 1rem;
}

.genome-metric {
    text-align: center;
    padding: 1rem;
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
}

.genome-metric.highlight {
    background: linear-gradient(135deg, rgba(129, 140, 248, 0.1), rgba(99, 102, 241, 0.1));
    border: 1px solid var(--accent-color);
}

.genome-metric-value {
    font-size: 1.75rem;
    font-weight: 700;
    color: var(--accent-color);
    line-height: 1.2;
}

.genome-metric-label {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-top: 0.25rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.genome-metric-hint {
    font-size: 0.75rem;
    color: var(--text-muted);
    margin-top: 0.25rem;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}

/* Genome Interpretation - Dark Theme */
.genome-interpretation {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem 1.25rem;
    margin-bottom: 1.5rem;
    border: 1px solid var(--border-color);
}

.genome-interpretation h4 {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.interp-items {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
}

.interp-item {
    display: flex;
    gap: 0.5rem;
    font-size: 0.8rem;
}

.interp-label {
    font-weight: 600;
    color: var(--text-primary);
    min-width: 160px;
    flex-shrink: 0;
}

.interp-desc {
    color: var(--text-secondary);
}

.plot-row {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
    gap: 1.5rem;
    margin-bottom: 1.5rem;
}

.plot-container {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-color);
}

.plot-container.full-width {
    grid-column: 1 / -1;
}

.plot-title {
    font-size: 0.95rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.table-container {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-color);
}

.table-controls {
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    margin-bottom: 1rem;
    padding-bottom: 1rem;
    border-bottom: 1px solid var(--border-color);
}

.table-search {
    flex: 1;
    min-width: 200px;
    padding: 0.5rem 0.75rem;
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    background: var(--bg-secondary);
    color: var(--text-primary);
}

.table-filter {
    padding: 0.5rem 0.75rem;
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    background: var(--bg-secondary);
    color: var(--text-primary);
}

.export-btn {
    padding: 0.5rem 1rem;
    background: var(--accent-color);
    color: white;
    border: none;
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    font-weight: 500;
    cursor: pointer;
}

.column-selector {
    position: relative;
}

.column-btn {
    padding: 0.5rem 1rem;
    background: var(--bg-tertiary);
    color: var(--text-primary);
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    font-size: 0.875rem;
    font-weight: 500;
    cursor: pointer;
    transition: background 0.2s ease;
}

.column-btn:hover {
    background: var(--bg-secondary);
}

.column-menu {
    display: none;
    position: absolute;
    top: 100%;
    right: 0;
    margin-top: 0.5rem;
    padding: 0.75rem;
    background: var(--bg-primary);
    border: 1px solid var(--border-color);
    border-radius: var(--radius-md);
    box-shadow: 0 4px 12px rgba(0,0,0,0.3);
    z-index: 100;
    min-width: 180px;
}

.column-menu.show {
    display: block;
}

.column-menu label {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.4rem 0.25rem;
    font-size: 0.85rem;
    cursor: pointer;
    white-space: nowrap;
    color: var(--text-primary);
}

.column-menu label:hover {
    background: var(--bg-secondary);
    border-radius: var(--radius-sm);
}

.column-menu input[type="checkbox"] {
    width: 16px;
    height: 16px;
    cursor: pointer;
}

.data-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.8rem;
}

.data-table th,
.data-table td {
    padding: 0.6rem 0.75rem;
    text-align: left;
    border-bottom: 1px solid var(--border-color);
}

.data-table th {
    background: var(--bg-secondary);
    font-weight: 600;
    color: var(--text-primary);
    cursor: pointer;
}

.data-table tbody tr:hover {
    background: var(--bg-tertiary);
}

.data-table .cell-known { color: var(--success-color); }
.data-table .cell-novel-species { color: var(--warning-color); }
.data-table .cell-novel-genus { color: var(--danger-color); }
.data-table .cell-conserved { color: var(--info-color); }

.table-pagination {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding-top: 1rem;
    border-top: 1px solid var(--border-color);
    margin-top: 1rem;
}

.page-info {
    font-size: 0.8rem;
    color: var(--text-secondary);
}

.page-btn {
    padding: 0.4rem 0.75rem;
    border: 1px solid var(--border-color);
    background: var(--bg-secondary);
    color: var(--text-primary);
    border-radius: var(--radius-sm);
    font-size: 0.8rem;
    cursor: pointer;
}

.page-number {
    font-size: 0.8rem;
    color: var(--text-secondary);
    padding: 0 0.5rem;
}

/* Data Tab Summary - Dark Theme */
.data-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1.25rem;
    margin-bottom: 1.25rem;
    border: 1px solid var(--border-color);
}

.data-intro h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 0.5rem 0;
}

.data-intro p {
    font-size: 0.875rem;
    color: var(--text-secondary);
    margin: 0;
    line-height: 1.5;
}

.data-stats {
    display: flex;
    gap: 1.5rem;
    margin-top: 1rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border-color);
}

.data-stat {
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 0.5rem 1rem;
    border-radius: var(--radius-sm);
    background: var(--bg-secondary);
}

.data-stat .stat-value {
    font-size: 1.25rem;
    font-weight: 700;
    color: var(--text-primary);
}

.data-stat .stat-label {
    font-size: 0.75rem;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.data-stat.known { border-left: 3px solid var(--success-color); }
.data-stat.novel { border-left: 3px solid var(--warning-color); }
.data-stat.uncertain { border-left: 3px solid var(--text-muted); }
.data-stat.off-target { border-left: 3px solid #3498db; }

/* Quick Filters - Dark Theme */
.quick-filters {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 1rem;
    padding: 0.75rem 1rem;
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
}

.filter-label {
    font-size: 0.8rem;
    color: var(--text-secondary);
    font-weight: 500;
    margin-right: 0.5rem;
}

.filter-chip {
    padding: 0.35rem 0.75rem;
    border: 1px solid var(--border-color);
    border-radius: 20px;
    background: var(--bg-primary);
    color: var(--text-primary);
    font-size: 0.75rem;
    cursor: pointer;
    transition: all 0.2s ease;
    white-space: nowrap;
}

.filter-chip:hover {
    border-color: var(--accent-color);
    background: var(--bg-tertiary);
}

.filter-chip.active {
    background: var(--accent-color);
    color: white;
    border-color: var(--accent-color);
}

.filter-chip.known.active { background: var(--success-color); border-color: var(--success-color); }
.filter-chip.novel-species.active { background: var(--warning-color); border-color: var(--warning-color); }
.filter-chip.novel-genus.active { background: var(--danger-color); border-color: var(--danger-color); }
.filter-chip.species-boundary.active { background: #8b5cf6; border-color: #8b5cf6; }
.filter-chip.ambiguous.active { background: var(--text-muted); border-color: var(--text-muted); }
.filter-chip.ambiguous-wg.active { background: #6b7280; border-color: #6b7280; }
.filter-chip.conserved.active { background: var(--info-color); border-color: var(--info-color); }
.filter-chip.unclassified.active { background: #d1d5db; border-color: #d1d5db; color: #1a1a2e; }
.filter-chip.off-target.active { background: #3498db; border-color: #3498db; }

/* Column Guide - Dark Theme */
.column-guide {
    margin-bottom: 1rem;
    background: var(--bg-primary);
    border: 1px solid var(--border-color);
    border-radius: var(--radius-md);
}

.column-guide summary {
    padding: 0.75rem 1rem;
    font-size: 0.85rem;
    font-weight: 500;
    color: var(--text-secondary);
    cursor: pointer;
    user-select: none;
}

.column-guide summary:hover {
    background: var(--bg-secondary);
}

.column-guide[open] summary {
    border-bottom: 1px solid var(--border-color);
}

.guide-content {
    padding: 0.75rem 1rem;
}

.guide-item {
    display: flex;
    gap: 1rem;
    padding: 0.5rem 0;
    font-size: 0.8rem;
    border-bottom: 1px solid var(--border-color);
}

.guide-item:last-child {
    border-bottom: none;
}

.guide-col {
    font-weight: 600;
    color: var(--text-primary);
    min-width: 100px;
}

.guide-desc {
    color: var(--text-secondary);
}

/* Enhanced Table Controls - Dark Theme */
.controls-row {
    display: flex;
    gap: 0.75rem;
    width: 100%;
    align-items: center;
}

.controls-right {
    display: flex;
    gap: 0.5rem;
    margin-left: auto;
}

.filter-status {
    display: flex;
    align-items: center;
    gap: 1rem;
    font-size: 0.85rem;
}

.results-count {
    color: var(--text-secondary);
}

.results-count strong {
    color: var(--text-primary);
    font-weight: 600;
}

.clear-btn {
    padding: 0.3rem 0.6rem;
    background: transparent;
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    font-size: 0.75rem;
    color: var(--text-secondary);
    cursor: pointer;
    transition: all 0.2s ease;
}

.clear-btn:hover {
    background: var(--danger-color);
    color: white;
    border-color: var(--danger-color);
}

.menu-header {
    font-size: 0.75rem;
    font-weight: 600;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.03em;
    padding: 0.25rem 0.25rem 0.5rem 0.25rem;
    border-bottom: 1px solid var(--border-color);
    margin-bottom: 0.5rem;
}

.table-wrapper {
    overflow-x: auto;
}

/* Methods Footer Panel */
.methods-footer {
    max-width: 1400px;
    margin: 1rem auto 0;
    padding: 0 1.25rem 1.25rem;
}

.report-footer {
    text-align: center;
    padding: 2rem;
    color: var(--text-muted);
    font-size: 0.8rem;
    border-top: 1px solid var(--border-color);
    margin-top: 2rem;
}

.report-footer a {
    color: var(--accent-color);
}

@media (max-width: 768px) {
    .plot-row { grid-template-columns: 1fr; }
    .metric-cards { grid-template-columns: repeat(2, 1fr); }
    .diversity-summary { padding: 1.25rem; }
    .diversity-cards { grid-template-columns: 1fr; gap: 1rem; }
    .diversity-pct { font-size: 2rem; }
    .breakdown-grid { grid-template-columns: 1fr; gap: 1rem; }
    .diversity-bar-legend { flex-wrap: wrap; gap: 1rem; }
    .dist-metrics { grid-template-columns: repeat(2, 1fr); }
    .interpretation-grid { grid-template-columns: 1fr; }
}
"""


# =============================================================================
# Utility Functions
# =============================================================================

# =============================================================================
# Novel Diversity Section Styles
# =============================================================================

NOVEL_SECTION_CSS: str = """
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
    background: var(--bg-primary);
    border: 1px solid var(--border-color);
    border-radius: 8px;
    padding: 1.25rem;
    text-align: center;
}

.novel-metric.highlight {
    border-color: var(--accent-color);
    background: linear-gradient(135deg, var(--bg-primary), rgba(102, 126, 234, 0.1));
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
    background: var(--bg-primary);
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
    background: var(--bg-secondary);
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
    background: var(--bg-primary);
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

/* Cluster table hint */
.cluster-table-hint {
    font-size: 0.8rem;
    color: var(--text-muted);
    margin-bottom: 0.75rem;
}

/* Neighborhood Panels */
.neighborhood-panel {
    padding: 1rem;
    background: var(--bg-secondary, #f8f9fa);
    border-top: 1px solid var(--border-color);
}

.neighborhood-panel-row td {
    padding: 0 !important;
    border-bottom: none !important;
}

.nbr-section { margin-bottom: 1rem; }
.nbr-section h5 {
    font-size: 0.8rem;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    margin-bottom: 0.5rem;
    color: var(--text-secondary);
}

.nbr-genera-table { width: 100%; font-size: 0.8rem; border-collapse: collapse; }
.nbr-genera-table th,
.nbr-genera-table td {
    padding: 0.3rem 0.5rem;
    text-align: left;
    border-bottom: 1px solid var(--border-color);
}
.nbr-genera-table th {
    font-weight: 600;
    color: var(--text-secondary);
}
.nbr-genera-table code {
    font-size: 0.75rem;
    background: var(--bg-tertiary, #e4e7eb);
    padding: 0.1rem 0.3rem;
    border-radius: 2px;
}

.nbr-metrics {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 0.5rem;
}

.nbr-metric {
    padding: 0.5rem;
    background: var(--bg-primary, white);
    border-radius: 4px;
    text-align: center;
    border: 1px solid var(--border-color);
}

.nbr-label {
    display: block;
    font-size: 0.7rem;
    color: var(--text-muted);
    text-transform: uppercase;
}

.nbr-value {
    display: block;
    font-size: 0.9rem;
    font-weight: 600;
}

/* Support score colors */
.support-tag {
    padding: 0.15rem 0.4rem;
    border-radius: 3px;
    font-size: 0.75rem;
    font-weight: 600;
}

.support-high { background: #dcfce7; color: #166534; }
.support-medium { background: #fef3c7; color: #92400e; }
.support-low { background: #fecaca; color: #991b1b; }

/* Context column */
.cluster-context {
    font-size: 0.8rem;
    max-width: 300px;
    color: var(--text-secondary);
}

@media (max-width: 768px) {
    .novel-metrics {
        grid-template-columns: repeat(2, 1fr);
    }

    .confidence-item {
        grid-template-columns: 1fr;
        gap: 0.5rem;
    }

    .phylo-item {
        grid-template-columns: 1fr;
        gap: 0.25rem;
    }

    .nbr-metrics {
        grid-template-columns: repeat(2, 1fr);
    }
}
"""

# =============================================================================
# Methods Section CSS
# =============================================================================

METHODS_SECTION_CSS: str = """
/* Methods Section Styles */
.methods-section {
    max-width: 900px;
    margin: 0 auto;
    line-height: 1.7;
}

.methods-intro {
    background: var(--card-bg);
    padding: 1.5rem;
    border-radius: 8px;
    margin-bottom: 2rem;
    border-left: 4px solid var(--primary-color);
}

.methods-intro p {
    margin: 0;
    color: var(--text-secondary);
}

.methods-intro a {
    color: var(--primary-color);
    text-decoration: none;
}

.methods-intro a:hover {
    text-decoration: underline;
}

.methods-toc {
    background: var(--card-bg);
    padding: 1.5rem;
    border-radius: 8px;
    margin-bottom: 2rem;
}

.methods-toc h3 {
    margin-top: 0;
    margin-bottom: 1rem;
    font-size: 1.1rem;
    color: var(--text-primary);
}

.methods-toc ol {
    margin: 0;
    padding-left: 1.5rem;
}

.methods-toc li {
    margin-bottom: 0.5rem;
}

.methods-toc a {
    color: var(--primary-color);
    text-decoration: none;
}

.methods-toc a:hover {
    text-decoration: underline;
}

.method-section {
    background: var(--card-bg);
    padding: 2rem;
    border-radius: 8px;
    margin-bottom: 1.5rem;
}

.method-section h3 {
    margin-top: 0;
    margin-bottom: 1.5rem;
    font-size: 1.3rem;
    color: var(--text-primary);
    border-bottom: 2px solid var(--primary-color);
    padding-bottom: 0.5rem;
}

.method-section h4 {
    margin-top: 1.5rem;
    margin-bottom: 1rem;
    font-size: 1.1rem;
    color: var(--text-primary);
}

.method-section p {
    color: var(--text-secondary);
    margin-bottom: 1rem;
}

.method-section ul, .method-section ol {
    color: var(--text-secondary);
    margin-bottom: 1rem;
    padding-left: 1.5rem;
}

.method-section li {
    margin-bottom: 0.5rem;
}

.method-highlight {
    background: rgba(var(--primary-rgb), 0.1);
    border-left: 4px solid var(--primary-color);
    padding: 1rem 1.5rem;
    margin: 1.5rem 0;
    border-radius: 0 8px 8px 0;
}

.method-highlight strong {
    color: var(--primary-color);
}

/* Formula boxes */
.formula-box {
    background: var(--bg-primary);
    border: 1px solid var(--border-color);
    border-radius: 8px;
    padding: 1.5rem;
    margin: 1.5rem 0;
    text-align: center;
}

.formula-box .formula-title {
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 1rem;
    text-align: left;
}

.formula-box .formula {
    font-family: 'Courier New', Consolas, monospace;
    font-size: 1.1rem;
    color: var(--text-primary);
    margin: 0.75rem 0;
    padding: 0.5rem;
    background: var(--card-bg);
    border-radius: 4px;
    display: inline-block;
}

.formula-box .formula-legend {
    font-size: 0.9rem;
    color: var(--text-muted);
    margin-top: 0.75rem;
    font-style: italic;
}

.formula-box .formula-table {
    width: 100%;
    margin-top: 1rem;
    text-align: left;
}

.formula-box .formula-table td {
    padding: 0.5rem;
    border-bottom: 1px solid var(--border-color);
    font-family: 'Courier New', Consolas, monospace;
    font-size: 0.95rem;
}

.formula-box .formula-table td:first-child {
    width: 30%;
    color: var(--text-muted);
}

.formula-box .formula-table td:last-child {
    color: var(--text-muted);
    font-family: inherit;
}

/* Methods tables */
.methods-table {
    width: 100%;
    border-collapse: collapse;
    margin: 1rem 0 1.5rem 0;
    font-size: 0.95rem;
}

.methods-table.small {
    max-width: 300px;
}

.methods-table th {
    background: var(--bg-primary);
    color: var(--text-primary);
    font-weight: 600;
    text-align: left;
    padding: 0.75rem 1rem;
    border-bottom: 2px solid var(--border-color);
}

.methods-table td {
    padding: 0.75rem 1rem;
    border-bottom: 1px solid var(--border-color);
    color: var(--text-secondary);
}

.methods-table tr:hover {
    background: var(--bg-hover);
}

.methods-table .row-known { background: rgba(34, 197, 94, 0.1); }
.methods-table .row-novel-species { background: rgba(234, 179, 8, 0.1); }
.methods-table .row-novel-genus { background: rgba(249, 115, 22, 0.1); }
.methods-table .row-unclassified { background: rgba(148, 163, 184, 0.1); }
.methods-table .row-confident { background: rgba(34, 197, 94, 0.1); }
.methods-table .row-boundary { background: rgba(234, 179, 8, 0.1); }
.methods-table .row-conserved { background: rgba(59, 130, 246, 0.1); }

.methods-table .priority-high { background: rgba(34, 197, 94, 0.15); }
.methods-table .priority-medium { background: rgba(234, 179, 8, 0.15); }
.methods-table .priority-low { background: rgba(249, 115, 22, 0.15); }
.methods-table .priority-uncertain { background: rgba(239, 68, 68, 0.15); }

.methods-note {
    font-size: 0.9rem;
    color: var(--text-muted);
    font-style: italic;
    margin-top: 1rem;
}

/* Decision tree visualization */
.decision-tree {
    font-family: 'Courier New', Consolas, monospace;
    padding: 1.5rem;
    background: var(--bg-primary);
    border-radius: 8px;
    margin: 1.5rem 0;
    line-height: 1.8;
}

.tree-node {
    padding: 0.5rem 1rem;
    background: var(--card-bg);
    border: 1px solid var(--border-color);
    border-radius: 4px;
    display: inline-block;
    margin: 0.25rem 0;
}

.tree-node.start {
    background: var(--primary-color);
    color: white;
    border-color: var(--primary-color);
}

.tree-node.result {
    font-weight: 600;
}

.tree-arrow {
    color: var(--text-muted);
    margin: 0.25rem 0;
    padding-left: 2rem;
}

.tree-branch {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin: 0.25rem 0;
}

.tree-condition {
    padding: 0.5rem 1rem;
    background: var(--bg-secondary);
    border: 1px dashed var(--border-color);
    border-radius: 4px;
}

.tree-yes {
    color: var(--text-muted);
}

/* Category badges in decision tree */
.cat-known { color: #22c55e; font-weight: 600; }
.cat-novel-species { color: #eab308; font-weight: 600; }
.cat-novel-genus { color: #f97316; font-weight: 600; }
.cat-boundary { color: #3b82f6; font-weight: 600; }
.cat-ambiguous { color: #6b7280; font-weight: 600; }
.cat-unclassified { color: #94a3b8; font-weight: 600; }

/* References list */
.references-list {
    padding-left: 1.5rem;
}

.references-list li {
    margin-bottom: 1rem;
    color: var(--text-secondary);
    line-height: 1.6;
}

.references-list strong {
    color: var(--text-primary);
}

.references-list em {
    color: var(--text-secondary);
}

.references-list a {
    color: var(--primary-color);
    text-decoration: none;
    font-size: 0.9rem;
}

.references-list a:hover {
    text-decoration: underline;
}

/* Code styling in methods */
.method-section code {
    font-family: 'Courier New', Consolas, monospace;
    background: var(--bg-primary);
    padding: 0.2rem 0.4rem;
    border-radius: 3px;
    font-size: 0.9rem;
    color: var(--primary-color);
}

/* Responsive adjustments */
@media (max-width: 768px) {
    .method-section {
        padding: 1.5rem;
    }

    .formula-box {
        padding: 1rem;
        overflow-x: auto;
    }

    .formula-box .formula {
        font-size: 0.95rem;
    }

    .decision-tree {
        font-size: 0.85rem;
        overflow-x: auto;
    }

    .tree-branch {
        flex-direction: column;
        align-items: flex-start;
    }

    .methods-table {
        font-size: 0.85rem;
    }

    .methods-table th,
    .methods-table td {
        padding: 0.5rem;
    }
}
"""


# =============================================================================
# Phylogeny Section Styles
# =============================================================================

PHYLOGENY_STYLES: str = """
/* Phylogeny Section */
.phylogeny-section {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-sm);
    border: 1px solid var(--border-light);
}

.phylogeny-header {
    margin-bottom: 1.5rem;
}

.phylogeny-header h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.5rem;
}

.phylogeny-description {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
}

/* Heatmap Toggle Controls (ANI/AAI representative filter) */
.heatmap-toggle-controls {
    display: flex;
    gap: 0;
    margin-bottom: 1rem;
}

.heatmap-toggle-btn {
    padding: 0.5rem 1.25rem;
    border: 1px solid var(--border-color);
    background: var(--bg-secondary);
    color: var(--text-secondary);
    font-size: 0.85rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s ease;
}

.heatmap-toggle-btn:first-child {
    border-radius: var(--radius-sm) 0 0 var(--radius-sm);
    border-right: none;
}

.heatmap-toggle-btn:last-child {
    border-radius: 0 var(--radius-sm) var(--radius-sm) 0;
}

.heatmap-toggle-btn:hover {
    background: var(--bg-tertiary);
}

.heatmap-toggle-btn.active {
    background: var(--accent-color);
    color: white;
    border-color: var(--accent-color);
}

/* Phylogeny Controls */
.phylogeny-controls {
    display: flex;
    gap: 0.75rem;
    margin-bottom: 1.5rem;
    flex-wrap: wrap;
}

.phylogeny-controls .btn {
    padding: 0.5rem 1rem;
    border: 1px solid var(--border-color);
    border-radius: var(--radius-sm);
    background: var(--bg-secondary);
    color: var(--text-primary);
    font-size: 0.85rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s ease;
}

.phylogeny-controls .btn:hover {
    background: var(--bg-tertiary);
    border-color: var(--accent-color);
}

.phylogeny-controls .btn:active {
    background: var(--accent-color);
    color: white;
}

.phylogeny-controls .btn-primary {
    background: var(--accent-color);
    color: white;
    border-color: var(--accent-color);
}

.phylogeny-controls .btn-primary:hover {
    background: var(--accent-hover);
}

/* Phylogeny Legend */
.phylogeny-legend {
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
    padding: 1rem;
    margin-bottom: 1.5rem;
}

.phylogeny-legend .legend-title {
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.phylogeny-legend .legend-items {
    display: flex;
    flex-wrap: wrap;
    gap: 1.5rem;
    margin-bottom: 0.75rem;
}

.phylogeny-legend .legend-item {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    font-size: 0.8rem;
    color: var(--text-secondary);
}

.phylogeny-legend .legend-circle {
    width: 14px;
    height: 14px;
    border-radius: 50%;
    border: 2px solid white;
    box-shadow: 0 1px 2px rgba(0,0,0,0.2);
}

.phylogeny-legend .legend-circle.reference {
    background-color: var(--success-color);
}

.phylogeny-legend .legend-circle.novel-species {
    background-color: var(--warning-color);
}

.phylogeny-legend .legend-circle.novel-genus {
    background-color: var(--danger-color);
}

.phylogeny-legend .legend-badges {
    display: flex;
    gap: 0.75rem;
    padding-top: 0.75rem;
    border-top: 1px solid var(--border-light);
}

/* Confidence Badges */
.confidence-badge {
    display: inline-block;
    padding: 0.2rem 0.5rem;
    border-radius: 10px;
    font-size: 0.7rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.confidence-badge.high {
    background-color: rgba(46, 204, 113, 0.2);
    color: var(--success-color);
    border: 1px solid var(--success-color);
}

.confidence-badge.medium {
    background-color: rgba(243, 156, 18, 0.2);
    color: var(--warning-color);
    border: 1px solid var(--warning-color);
}

.confidence-badge.low {
    background-color: rgba(231, 76, 60, 0.2);
    color: var(--danger-color);
    border: 1px solid var(--danger-color);
}

/* Tree Container */
.phylogeny-tree-container {
    background: var(--bg-secondary);
    border-radius: var(--radius-md);
    padding: 1rem;
    min-height: 400px;
    overflow: auto;
    border: 1px solid var(--border-light);
}

.phylogeny-tree-container svg {
    display: block;
    margin: 0 auto;
}

/* Tree Links */
.phylogeny-tree-container .link {
    fill: none;
    stroke: var(--border-color);
    stroke-width: 1.5px;
}

.phylogeny-tree-container .link.novel {
    stroke: var(--warning-color);
    stroke-dasharray: 5, 3;
    stroke-width: 2px;
}

/* Tree Nodes */
.phylogeny-tree-container .node circle {
    stroke: white;
    stroke-width: 1.5px;
}

.phylogeny-tree-container .node.leaf text {
    font-size: 10px;
    fill: var(--text-primary);
}

.phylogeny-tree-container .node.novel circle {
    stroke: var(--warning-color);
    stroke-width: 2px;
}

/* Phylogeny Tooltip */
.phylogeny-tooltip {
    position: absolute;
    background: var(--bg-primary);
    border: 1px solid var(--border-color);
    border-radius: var(--radius-md);
    padding: 0.75rem 1rem;
    box-shadow: var(--shadow-md);
    font-size: 0.8rem;
    color: var(--text-primary);
    max-width: 280px;
    z-index: 1000;
    pointer-events: none;
}

.phylogeny-tooltip strong {
    display: block;
    margin-bottom: 0.35rem;
    color: var(--text-primary);
}

.phylogeny-tooltip br {
    line-height: 1.8;
}

/* Responsive adjustments */
@media (max-width: 768px) {
    .phylogeny-controls {
        flex-direction: column;
    }

    .phylogeny-controls .btn {
        width: 100%;
        text-align: center;
    }

    .phylogeny-legend .legend-items {
        flex-direction: column;
        gap: 0.5rem;
    }

    .phylogeny-tree-container {
        min-height: 300px;
    }
}
"""


BORDERLINE_ANALYSIS_CSS: str = """
.borderline-cards {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
    margin: 1.5rem 0;
}

.borderline-card {
    background: var(--bg-secondary);
    border: 1px solid var(--border-color);
    border-radius: var(--radius-md);
    padding: 1rem;
    text-align: center;
}

.borderline-card .card-count {
    font-size: 1.8rem;
    font-weight: 700;
    color: var(--text-primary);
    line-height: 1.2;
}

.borderline-card .card-pct {
    font-size: 0.85rem;
    color: var(--text-secondary);
    margin-bottom: 0.3rem;
}

.borderline-card .card-label {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--text-secondary);
    margin-bottom: 0.25rem;
}

.borderline-card .card-desc {
    font-size: 0.7rem;
    color: var(--text-muted);
}

.borderline-card.highlight-amber {
    border-color: var(--warning-color);
    background: rgba(243, 156, 18, 0.06);
}

.borderline-card.highlight-red {
    border-color: var(--danger-color);
    background: rgba(231, 76, 60, 0.06);
}
"""


def get_css_styles(theme: str = "light") -> str:
    """
    Get CSS styles for specified theme.

    Args:
        theme: Theme name ('light' or 'dark')

    Returns:
        CSS stylesheet string
    """
    base_theme = DARK_THEME if theme.lower() == "dark" else LIGHT_THEME
    return base_theme + NOVEL_SECTION_CSS + METHODS_SECTION_CSS + PHYLOGENY_STYLES + BORDERLINE_ANALYSIS_CSS
