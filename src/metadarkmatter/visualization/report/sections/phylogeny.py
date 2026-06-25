"""ReportGenerator phylogeny section."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

from metadarkmatter.visualization.report.sections._base import _ReportBase

if TYPE_CHECKING:
    import pandas as pd



logger = logging.getLogger(__name__)


class PhylogenyMixin(_ReportBase):
    """Phylogeny report-section builders."""

    def _build_phylogeny_content(self) -> str | None:
        """Build phylogeny HTML content without tab wrapper.

        Generates an interactive phylogenetic tree from the ANI matrix (or
        a user-provided tree) with novel cluster placement.

        Returns:
            HTML string with phylogeny visualization, or None if the tree
            cannot be built (e.g., skip requested, no ANI matrix, or fewer
            than 3 genomes).
        """
        if self.config.skip_phylogeny:
            return None
        if self.ani_matrix is None or len(self.ani_matrix) < 3:
            return None

        ani_pd = self.ani_matrix.to_pandas()
        if "genome" in ani_pd.columns:
            ani_pd = ani_pd.set_index("genome")

        return self._build_phylogeny_section(
            ani_pd,
            user_tree_path=self.config.user_tree_path,
        )

    def _build_phylogeny_section(
        self,
        ani_matrix_pd: pd.DataFrame | None,
        user_tree_path: Path | None = None,
    ) -> str | None:
        """Build interactive phylogenetic tree section.

        Generates a phylogenetic tree from the ANI matrix (or loads a user-provided
        tree) and places novel clusters at estimated positions. The tree is rendered
        using D3.js for interactive exploration.

        Args:
            ani_matrix_pd: Pandas DataFrame with ANI values. Rows and columns
                should be genome accessions. Should be pre-indexed (no 'genome'
                column as index).
            user_tree_path: Optional path to user-provided Newick tree file.
                If provided, this tree is used instead of building from ANI.

        Returns:
            HTML string with phylogeny section content and D3.js visualization
            code, or None if the tree cannot be built (e.g., fewer than 3 genomes).
        """
        import json

        from metadarkmatter.visualization.report.templates import (
            PHYLOGENY_SECTION_TEMPLATE,
            PHYLOTREE_JS_TEMPLATE,
        )

        if ani_matrix_pd is None or len(ani_matrix_pd) < 3:
            return None

        try:
            from metadarkmatter.core.phylogeny import (
                ani_to_newick,
                extract_novel_clusters,
                load_user_tree,
                place_novel_clusters,
            )

            # Build or load tree
            newick: str | None
            if user_tree_path is not None:
                newick = load_user_tree(user_tree_path, ani_matrix_pd)
                tree_source_note = "Tree source: user-provided Newick file."
            else:
                newick = ani_to_newick(ani_matrix_pd)
                tree_source_note = "Tree source: neighbor-joining from ANI matrix."

            if newick is None:
                return None

            # Extract novel clusters from classifications
            novel_clusters = extract_novel_clusters(self.df, min_reads=3)

            # Place novel clusters on the tree
            if novel_clusters:
                final_newick = place_novel_clusters(
                    newick, novel_clusters, ani_matrix_pd
                )
            else:
                final_newick = newick

            # Strip BioPython comment blocks from Newick (JavaScript parser can't handle them)
            # Comments look like: NSP_001:0.1[{"classification": "novel_species", ...}]
            # We already have this data in the annotations dict
            import re
            final_newick = re.sub(r'\[[^\]]*\]', '', final_newick)

            # Build annotations dict for novel clusters
            annotations: dict[str, dict[str, Any]] = {}
            for cluster in novel_clusters:
                # Novel clusters are named by their cluster_id
                node_name = cluster.cluster_id
                annotations[node_name] = {
                    "type": cluster.classification,
                    "read_count": cluster.read_count,
                    "mean_novelty": cluster.mean_novelty,
                    "mean_uncertainty": cluster.mean_uncertainty,
                    "confidence": cluster.confidence_rating,
                    "nearest_ref": cluster.best_match_genome,
                    "est_ani": cluster.mean_identity,
                    "is_novel": True,
                }

            # Build genome labels for tree leaf display
            genome_labels = {
                acc: self._get_genome_label(acc, max_species_len=20)
                for acc in ani_matrix_pd.columns
            }

            # Build tree data JSON for D3.js visualization
            tree_data = {
                "newick": final_newick,
                "annotations": annotations,
                "genome_labels": genome_labels,
                "tip_count": len(ani_matrix_pd) + len(novel_clusters),
                "novel_count": len(novel_clusters),
                "source_note": tree_source_note,
            }

            section_html = PHYLOGENY_SECTION_TEMPLATE.format(
                tree_source_note=tree_source_note,
                tree_data_json=json.dumps(tree_data),
            )

            # PHYLOTREE_JS_TEMPLATE uses {{ and }} for JS braces, call .format() to convert
            from metadarkmatter.visualization.report.assets import get_d3_script_tag
            d3_script_tag = get_d3_script_tag(self.config.report_mode)
            return section_html + PHYLOTREE_JS_TEMPLATE.format(d3_script_tag=d3_script_tag)

        except ImportError as e:
            logger.warning(f"Phylogeny dependencies not available: {e}")
            return None
        except Exception as e:
            logger.warning(f"Error building phylogeny section: {e}")
            return None
