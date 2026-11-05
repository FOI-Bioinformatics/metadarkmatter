"""Shared builders for the HTML report unit tests."""

from __future__ import annotations

import json
import re
from pathlib import Path

import polars as pl

GENOMES = ["GCF_000001.1", "GCF_000002.1", "GCF_000003.1", "GCF_000004.1"]

CALLS = [
    "Known Species",
    "Novel Species",
    "Novel Genus",
    "Species Boundary",
    "Ambiguous",
    "Ambiguous Within Genus",
    "Conserved Region",
    "Unclassified",
    "Off-target",
]


def make_classifications(n: int = 90) -> pl.DataFrame:
    """Deterministic classification frame covering every category.

    Includes ``legacy_call``, ``uncertainty_type`` and ``confidence_score`` so
    the call-comparison / run-info builders have data to work with.
    """
    calls = [CALLS[i % len(CALLS)] for i in range(n)]
    legacy = [calls[(i + (1 if i % 7 == 0 else 0)) % n] for i in range(n)]
    genomes = [GENOMES[i % len(GENOMES)] for i in range(n)]
    identity = [round(100.0 - (i % 30) * 0.9, 2) for i in range(n)]
    unc: list[float | None] = [
        None if c == "Off-target" else round((i % 11) * 0.7, 2) for i, c in enumerate(calls)
    ]
    return pl.DataFrame(
        {
            "read_id": [f"read_{i:04d}/1" for i in range(n)],
            "best_match_genome": genomes,
            "top_hit_identity": identity,
            "novelty_index": [round(100.0 - x, 2) for x in identity],
            "placement_uncertainty": unc,
            "num_ambiguous_hits": [1 + (i % 3) for i in range(n)],
            "identity_gap": [None if i % 4 == 0 else round(i % 9 * 0.5, 2) for i in range(n)],
            "taxonomic_call": calls,
            "legacy_call": legacy,
            "is_novel": [c in ("Novel Species", "Novel Genus") for c in calls],
            "uncertainty_type": ["measured" if i % 3 else "inferred" for i in range(n)],
            "confidence_score": [round(50 + (i % 40), 1) for i in range(n)],
            "posterior_entropy": [round((i % 20) * 0.12, 3) for i in range(n)],
            "p_known_species": [0.5] * n,
            "p_novel_species": [0.2] * n,
            "p_novel_genus": [0.1] * n,
            "p_species_boundary": [0.1] * n,
            "p_ambiguous": [0.05] * n,
            "p_unclassified": [0.05] * n,
        },
        schema_overrides={"placement_uncertainty": pl.Float64},
    )


def make_ani() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "genome": GENOMES,
            "GCF_000001.1": [100.0, 95.5, 85.2, 78.0],
            "GCF_000002.1": [95.5, 100.0, 84.1, 77.5],
            "GCF_000003.1": [85.2, 84.1, 100.0, 92.0],
            "GCF_000004.1": [78.0, 77.5, 92.0, 100.0],
        }
    )


def make_metadata():
    from metadarkmatter.core.metadata import GenomeMetadata

    return GenomeMetadata(
        pl.DataFrame(
            {
                "accession": GENOMES,
                "species": [
                    "Francisella tularensis",
                    "Francisella novicida",
                    "Allofrancisella sp",
                    "Allofrancisella sp",
                ],
                "genus": ["Francisella", "Francisella", "Allofrancisella", "Allofrancisella"],
                "family": ["Francisellaceae"] * 4,
                "gtdb_taxonomy": ["d__Bacteria"] * 4,
            }
        )
    )


def render(tmp_path: Path, df: pl.DataFrame | None = None, **kwargs) -> str:
    """Render a report with the default fixtures; returns the HTML text."""
    from metadarkmatter.visualization.report.generator import ReportGenerator

    config = kwargs.pop("config", None)
    gen = ReportGenerator(
        df if df is not None else make_classifications(),
        config=config,
        ani_matrix=kwargs.pop("ani_matrix", make_ani()),
        genome_metadata=kwargs.pop("genome_metadata", make_metadata()),
        **kwargs,
    )
    out = tmp_path / "report.html"
    gen.generate(out)
    return out.read_text(encoding="utf-8")


def embedded_json(html_text: str, element_id: str) -> object:
    """Parse the JSON of ``<script type="application/json" id=element_id>``."""
    m = re.search(
        r'<script type="application/json" id="' + re.escape(element_id) + r'">(.*?)</script>',
        html_text,
        re.S,
    )
    assert m is not None, f"no embedded JSON with id {element_id!r}"
    return json.loads(m.group(1))
