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
    --bg-primary: #ffffff;
    --bg-secondary: #f8f9fa;
    --bg-tertiary: #e9ecef;
    --text-primary: #212529;
    --text-secondary: #6c757d;
    --text-muted: #adb5bd;
    --border-color: #dee2e6;
    --border-light: #e9ecef;
    --accent-color: #0d6efd;
    --accent-hover: #0b5ed7;
    --success-color: #2ecc71;
    --warning-color: #f39c12;
    --danger-color: #e74c3c;
    --info-color: #95a5a6;
    --shadow-sm: 0 1px 2px rgba(0,0,0,0.05);
    --shadow-md: 0 4px 6px rgba(0,0,0,0.07);
    --shadow-lg: 0 10px 15px rgba(0,0,0,0.1);
    --radius-sm: 4px;
    --radius-md: 8px;
    --radius-lg: 12px;
}

* {
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    background-color: var(--bg-secondary);
    color: var(--text-primary);
    line-height: 1.6;
    font-size: 14px;
}

/* Header */
.report-header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 2rem 2rem 1.5rem;
    text-align: center;
}

.report-header h1 {
    font-size: 1.75rem;
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
    gap: 0.25rem;
    background: var(--bg-primary);
    border-bottom: 1px solid var(--border-color);
    padding: 0.75rem 1rem;
    position: sticky;
    top: 0;
    z-index: 100;
    box-shadow: var(--shadow-sm);
}

.tab-btn {
    padding: 0.5rem 1.25rem;
    border: none;
    background: transparent;
    cursor: pointer;
    font-size: 0.875rem;
    font-weight: 500;
    color: var(--text-secondary);
    border-radius: var(--radius-sm);
    transition: all 0.2s ease;
}

.tab-btn:hover {
    background: var(--bg-secondary);
    color: var(--text-primary);
}

.tab-btn.active {
    background: var(--accent-color);
    color: white;
}

/* Main Content */
.report-content {
    max-width: 1400px;
    margin: 0 auto;
    padding: 1.5rem;
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
    font-size: 1.25rem;
    font-weight: 600;
    margin-bottom: 1.5rem;
    color: var(--text-primary);
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--accent-color);
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

/* Diversity Summary (Hero Section) */
.diversity-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 2rem;
    margin-bottom: 2rem;
    box-shadow: var(--shadow-md);
    border: 1px solid var(--border-light);
}

.diversity-headline {
    text-align: center;
    margin-bottom: 2rem;
}

.headline-value {
    font-size: 4rem;
    font-weight: 800;
    line-height: 1;
    margin-bottom: 0.5rem;
}

.headline-value.novel { color: var(--danger-color); }
.headline-value.known { color: var(--success-color); }

.headline-label {
    font-size: 1.5rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.25rem;
}

.headline-detail {
    font-size: 0.95rem;
    color: var(--text-secondary);
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

.diversity-pct {
    font-size: 2.5rem;
    font-weight: 700;
    line-height: 1.1;
}

.diversity-card.known .diversity-pct { color: var(--success-color); }
.diversity-card.novel .diversity-pct { color: var(--danger-color); }
.diversity-card.uncertain .diversity-pct { color: var(--info-color); }

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
    height: 24px;
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

.diversity-bar-legend {
    display: flex;
    justify-content: center;
    gap: 2rem;
    font-size: 0.85rem;
    color: var(--text-secondary);
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

/* Genome Highlights */
.genome-highlights {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem 1.25rem;
    margin-bottom: 1.5rem;
    border: 1px solid var(--border-light);
}

.genome-highlights h4 {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.highlights-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1rem;
}

.highlight-card {
    display: flex;
    align-items: flex-start;
    gap: 0.75rem;
    padding: 0.75rem;
    border-radius: var(--radius-md);
    background: var(--bg-secondary);
}

.highlight-icon {
    width: 32px;
    height: 32px;
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    font-weight: 700;
    font-size: 0.9rem;
    flex-shrink: 0;
}

.highlight-card.top .highlight-icon {
    background: var(--accent-color);
    color: white;
}

.highlight-card.novel .highlight-icon {
    background: var(--danger-color);
    color: white;
}

.highlight-card.confident .highlight-icon {
    background: var(--success-color);
    color: white;
}

.highlight-content {
    flex: 1;
    min-width: 0;
}

.highlight-title {
    font-size: 0.7rem;
    font-weight: 600;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.highlight-species {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-top: 0.25rem;
    font-style: italic;
}

.highlight-accession {
    font-size: 0.75rem;
    font-family: monospace;
    color: var(--text-muted);
    margin-top: 0.1rem;
}

.highlight-stats {
    display: flex;
    gap: 0.75rem;
    margin-top: 0.35rem;
}

.highlight-stats .stat {
    font-size: 0.75rem;
    color: var(--text-secondary);
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
.filter-chip.ambiguous.active { background: var(--text-muted); border-color: var(--text-muted); }
.filter-chip.conserved.active { background: var(--info-color); border-color: var(--info-color); }

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

    .headline-value {
        font-size: 2.5rem;
    }

    .headline-label {
        font-size: 1.1rem;
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

    .highlights-grid {
        grid-template-columns: 1fr;
    }

    .interp-item {
        flex-direction: column;
        gap: 0.25rem;
    }

    .interp-label {
        min-width: auto;
    }
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
        background: #667eea;
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
    --bg-primary: #1a1a2e;
    --bg-secondary: #16213e;
    --bg-tertiary: #0f3460;
    --text-primary: #eaeaea;
    --text-secondary: #94a3b8;
    --text-muted: #64748b;
    --border-color: #334155;
    --border-light: #1e293b;
    --accent-color: #818cf8;
    --accent-hover: #6366f1;
    --success-color: #34d399;
    --warning-color: #fbbf24;
    --danger-color: #f87171;
    --info-color: #94a3b8;
    --shadow-sm: 0 1px 2px rgba(0,0,0,0.3);
    --shadow-md: 0 4px 6px rgba(0,0,0,0.4);
    --shadow-lg: 0 10px 15px rgba(0,0,0,0.5);
    --radius-sm: 4px;
    --radius-md: 8px;
    --radius-lg: 12px;
}

* {
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    background-color: var(--bg-secondary);
    color: var(--text-primary);
    line-height: 1.6;
    font-size: 14px;
}

.report-header {
    background: linear-gradient(135deg, #4f46e5 0%, #7c3aed 100%);
    color: white;
    padding: 2rem 2rem 1.5rem;
    text-align: center;
}

.report-header h1 {
    font-size: 1.75rem;
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
    gap: 0.25rem;
    background: var(--bg-primary);
    border-bottom: 1px solid var(--border-color);
    padding: 0.75rem 1rem;
    position: sticky;
    top: 0;
    z-index: 100;
}

.tab-btn {
    padding: 0.5rem 1.25rem;
    border: none;
    background: transparent;
    cursor: pointer;
    font-size: 0.875rem;
    font-weight: 500;
    color: var(--text-secondary);
    border-radius: var(--radius-sm);
    transition: all 0.2s ease;
}

.tab-btn:hover {
    background: var(--bg-tertiary);
    color: var(--text-primary);
}

.tab-btn.active {
    background: var(--accent-color);
    color: white;
}

.report-content {
    max-width: 1400px;
    margin: 0 auto;
    padding: 1.5rem;
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
    font-size: 1.25rem;
    font-weight: 600;
    margin-bottom: 1.5rem;
    color: var(--text-primary);
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--accent-color);
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

/* Diversity Summary (Hero Section) - Dark Theme */
.diversity-summary {
    background: var(--bg-primary);
    border-radius: var(--radius-lg);
    padding: 2rem;
    margin-bottom: 2rem;
    box-shadow: var(--shadow-md);
    border: 1px solid var(--border-color);
}

.diversity-headline {
    text-align: center;
    margin-bottom: 2rem;
}

.headline-value {
    font-size: 4rem;
    font-weight: 800;
    line-height: 1;
    margin-bottom: 0.5rem;
}

.headline-value.novel { color: var(--danger-color); }
.headline-value.known { color: var(--success-color); }

.headline-label {
    font-size: 1.5rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.25rem;
}

.headline-detail {
    font-size: 0.95rem;
    color: var(--text-secondary);
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

.diversity-pct {
    font-size: 2.5rem;
    font-weight: 700;
    line-height: 1.1;
}

.diversity-card.known .diversity-pct { color: var(--success-color); }
.diversity-card.novel .diversity-pct { color: var(--danger-color); }
.diversity-card.uncertain .diversity-pct { color: var(--info-color); }

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
    height: 24px;
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

.diversity-bar-legend {
    display: flex;
    justify-content: center;
    gap: 2rem;
    font-size: 0.85rem;
    color: var(--text-secondary);
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

/* Genome Highlights - Dark Theme */
.genome-highlights {
    background: var(--bg-primary);
    border-radius: var(--radius-md);
    padding: 1rem 1.25rem;
    margin-bottom: 1.5rem;
    border: 1px solid var(--border-color);
}

.genome-highlights h4 {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.75rem;
}

.highlights-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 1rem;
}

.highlight-card {
    display: flex;
    align-items: flex-start;
    gap: 0.75rem;
    padding: 0.75rem;
    border-radius: var(--radius-md);
    background: var(--bg-secondary);
}

.highlight-icon {
    width: 32px;
    height: 32px;
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    font-weight: 700;
    font-size: 0.9rem;
    flex-shrink: 0;
}

.highlight-card.top .highlight-icon {
    background: var(--accent-color);
    color: white;
}

.highlight-card.novel .highlight-icon {
    background: var(--danger-color);
    color: white;
}

.highlight-card.confident .highlight-icon {
    background: var(--success-color);
    color: white;
}

.highlight-content {
    flex: 1;
    min-width: 0;
}

.highlight-title {
    font-size: 0.7rem;
    font-weight: 600;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.03em;
}

.highlight-species {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-top: 0.25rem;
    font-style: italic;
}

.highlight-accession {
    font-size: 0.75rem;
    font-family: monospace;
    color: var(--text-muted);
    margin-top: 0.1rem;
}

.highlight-stats {
    display: flex;
    gap: 0.75rem;
    margin-top: 0.35rem;
}

.highlight-stats .stat {
    font-size: 0.75rem;
    color: var(--text-secondary);
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
.filter-chip.ambiguous.active { background: var(--text-muted); border-color: var(--text-muted); }
.filter-chip.conserved.active { background: var(--info-color); border-color: var(--info-color); }

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
    .headline-value { font-size: 2.5rem; }
    .headline-label { font-size: 1.1rem; }
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

def get_css_styles(theme: str = "light") -> str:
    """
    Get CSS styles for specified theme.

    Args:
        theme: Theme name ('light' or 'dark')

    Returns:
        CSS stylesheet string
    """
    if theme.lower() == "dark":
        return DARK_THEME
    return LIGHT_THEME
