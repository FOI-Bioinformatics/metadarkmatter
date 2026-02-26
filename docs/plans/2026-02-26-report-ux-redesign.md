# Report UX Redesign: Dashboard Consolidation

**Date:** 2026-02-26
**Status:** Approved

## Problem

The current HTML report has 10 tabs with scattered information, poor visual hierarchy, and a dense layout that makes it hard to quickly answer "is there novel diversity?" The target audience is expert bioinformaticians who need fast signal extraction.

## Design Goals

- Answer "is there novel diversity?" in under 2 seconds
- Reduce tabs from 10 to 4-5
- Dense, Grafana-style dashboard aesthetic
- Clear visual hierarchy with compact KPI strips
- Collapsible panels for deep-dive content

## Tab Structure

```
[Summary] [Classification] [Novel Diversity*] [Reference] [Data]
```

*Novel Diversity only shown when novel reads exist.
Methods become a collapsible footer panel below all tabs.

### Merged Content Mapping

| Current Tab | New Location |
|-------------|-------------|
| Overview | Summary |
| Distributions | Classification |
| Species & Genomes | Reference (species bar) + Summary (top species mini-table) |
| Novel Diversity | Novel Diversity (unchanged) |
| Family Validation | Summary (key findings alerts + Off-target KPI) |
| Reference ANI | Reference (side-by-side with AAI) |
| Reference AAI | Reference (side-by-side with ANI) |
| Phylogeny | Reference (collapsible panel) |
| Classification Confidence | Classification (collapsible panel) |
| Recruitment | Reference (collapsible panel, if data provided) |
| Data | Data (unchanged) |
| Methods | Footer panel (collapsible) |

---

## Tab Designs

### 1. Summary Tab

Goal: answer "is there novel diversity?" immediately.

```
+------------------------------------------------------------------+
| KPI STRIP (horizontal, full-width)                                |
| [Total Reads] [Known %] [Novel %] [Uncertain %] [Off-target %*] |
+------------------------------------------------------------------+
|                                                                    |
| DIVERSITY BAR (full-width stacked bar, prominent)                 |
| [===Known===|==Novel==|=Uncertain=|Off-target]                    |
|                                                                    |
+-------------------------------+----------------------------------+
| KEY FINDINGS (left 60%)       | TOP SPECIES (right 40%)          |
| - Alert cards with icons      | Mini table: species, reads, %    |
| - "67% off-target detected"   | Top 8 species, sorted by count   |
| - "3 high-conf novel clusters"|                                  |
+-------------------------------+----------------------------------+
| CATEGORY BREAKDOWN (compact 2-row grid)                           |
| Known Species | Novel Species | Novel Genus | Species Boundary   |
| Ambiguous     | Conserved     | Unclassified| Off-target         |
+------------------------------------------------------------------+
```

Key changes from current Overview:
- KPI strip replaces the large hero headline (saves vertical space)
- Diversity bar is wider and more prominent
- Key Findings + Top Species side-by-side
- Category breakdown as a compact grid
- Sunburst moved to Classification tab
- Family Validation absorbed into KPI and Key Findings

### 2. Classification Tab

Goal: "how were reads classified and how confident are we?"

```
+------------------------------------------------------------------+
| METRIC STRIP                                                      |
| [Mean Novelty] [Mean Uncertainty] [Mean Entropy] [High-Conf %]   |
+------------------------------------------------------------------+
+-------------------------------+----------------------------------+
| NOVELTY vs UNCERTAINTY        | SUNBURST                         |
| SCATTER (left 60%)            | (right 40%)                      |
| - 2D scatter with regions     | - Diversity -> Category           |
| - Color by taxonomic_call     | - Interactive drill-down         |
+-------------------------------+----------------------------------+
+-------------------------------+----------------------------------+
| NOVELTY HISTOGRAM (left 50%)  | UNCERTAINTY HISTOGRAM (right 50%)|
+-------------------------------+----------------------------------+
+------------------------------------------------------------------+
| CONFIDENCE ANALYSIS (collapsible, default collapsed)              |
| Posterior entropy dist | Uncertainty type breakdown               |
| Borderline analysis    | MAP vs Legacy agreement                  |
+------------------------------------------------------------------+
```

Key changes:
- Scatter plot (primary diagnostic) is prominent and large
- Sunburst moves here from Overview
- Histograms side-by-side
- Confidence/Bayesian as collapsible deep-dive
- Confidence vs Novelty scatter dropped (redundant)

### 3. Novel Diversity Tab (conditional)

Goal: "what novel organisms did we find?"

```
+------------------------------------------------------------------+
| METRIC STRIP                                                      |
| [Clusters] [Novel Species] [Novel Genus] [High Confidence]       |
+------------------------------------------------------------------+
+------------------------------------------------------------------+
| CLUSTER TABLE (full-width, hero element)                          |
| ID | Type | Reads | Placement | Est. ANI | Uncertainty | Conf    |
+------------------------------------------------------------------+
+-------------------------------+----------------------------------+
| CLUSTER SCATTER (left 50%)    | PHYLOGENETIC SUNBURST (right 50%)|
+-------------------------------+----------------------------------+
+------------------------------------------------------------------+
| PHYLOGENETIC HEATMAP (collapsible)                                |
| Extended ANI/AAI matrix with novel clusters among references      |
+------------------------------------------------------------------+
```

Key changes:
- Cluster table is the hero element (most actionable)
- Confidence guide text removed
- Heatmap collapsible

### 4. Reference Tab

Goal: "what's the database and species landscape?"

```
+------------------------------------------------------------------+
| SPECIES BAR CHART (full-width, compact)                           |
+------------------------------------------------------------------+
+-------------------------------+----------------------------------+
| ANI HEATMAP (left 50%)        | AAI HEATMAP (right 50%)          |
| or full-width if only ANI     | (only if AAI matrix provided)    |
+-------------------------------+----------------------------------+
+------------------------------------------------------------------+
| PHYLOGENETIC TREE (collapsible)                                   |
| Interactive D3.js tree with novel cluster placement               |
+------------------------------------------------------------------+
+------------------------------------------------------------------+
| RECRUITMENT PLOTS (collapsible, only if data provided)            |
+------------------------------------------------------------------+
```

Key changes:
- Species bar chart moves from its own tab
- ANI + AAI side-by-side
- Phylogeny + Recruitment in collapsible panels
- Genome interpretation text removed

### 5. Data Tab

Unchanged. The interactive table with filters, search, pagination, and export is well-built.

### 6. Methods (Footer Panel)

Collapsible panel below all tab content. Default collapsed. Contains methodology documentation.

---

## Visual Style

### Color System

```css
--bg-page:       #f0f2f5    /* Cool gray page background */
--bg-card:       #ffffff    /* White cards */
--border-color:  #e4e7eb    /* Subtle card borders */
--header-bg:     #1a1a2e    /* Dark navy header */
--accent:        #667eea    /* Purple-blue accent */
--text-primary:  #1a1a2e
--text-secondary: #6b7280
--text-muted:    #9ca3af
```

Category colors (unchanged):
- Known: #22c55e (green)
- Novel Species: #f59e0b (amber)
- Novel Genus: #ef4444 (red)
- Uncertain: #94a3b8 (slate)
- Off-target: #3498db (blue)

### KPI Strip

Each metric card:
- Value: 1.75rem, weight 700, tabular-nums
- Label: 0.75rem, uppercase, letter-spacing 0.5px
- Left border: 3px, category color
- White card on gray page background

### Tab Navigation

Underline-style tab bar (replaces pill buttons):
- Active: bold + 2px bottom border in accent color
- Inactive: text-secondary
- Font: 0.8rem, weight 500

### Typography

- Numbers: `font-variant-numeric: tabular-nums`
- Section headers: 0.85rem, uppercase, letter-spacing 1px, weight 600
- Body: 0.8125rem (13px), line-height 1.5
- Tab buttons: 0.8rem, weight 500

### Spacing

- Page padding: 1.25rem
- Card gap: 0.75rem
- Section spacing: 1.5rem
- Card padding: 1rem

### Collapsible Panels

- Default collapsed with chevron + title
- Click to expand with smooth animation (max-height transition)
- Background: #f8f9fa for distinction from cards

---

## What Stays Unchanged

- Plotly chart library and chart types (scatter, histogram, heatmap, sunburst)
- Data table implementation (filters, sort, pagination, export)
- Plotly CDN loading strategy
- Column visibility localStorage persistence
- Responsive breakpoints (mobile at 768px)
- Conditional tab display logic
- Chart subsampling for performance

## What Gets Removed

- Hero headline (replaced by KPI strip)
- Confidence vs Novelty scatter (redundant)
- Genome highlights/interpretation cards (text-heavy, low value)
- Confidence guide text in Novel Diversity (obvious to experts)
- Methods as a tab (becomes footer panel)
- Purple gradient header (replaced by dark navy)
- Pill-style tab buttons (replaced by underline tabs)
