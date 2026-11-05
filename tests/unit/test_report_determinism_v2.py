"""Cross-process determinism of the new report sections (hash-seed independence)."""

from __future__ import annotations

import os
import re
import subprocess
import sys
import textwrap

import pytest

pytest.importorskip("plotly")

_SCRIPT = textwrap.dedent('''
    import sys
    from pathlib import Path
    from metadarkmatter.visualization.report.generator import ReportGenerator
    from metadarkmatter.visualization.report.models import ReportConfig
    from metadarkmatter.visualization.report.provenance import load_provenance
    from tests.unit.report_helpers import make_ani, make_classifications, make_metadata

    prov = load_provenance(Path(sys.argv[2]))
    ReportGenerator(
        make_classifications(300),
        config=ReportConfig(max_table_rows=120, input_names={"classifications": "c.csv"}),
        ani_matrix=make_ani(),
        genome_metadata=make_metadata(),
        provenance=prov,
    ).generate(Path(sys.argv[1]))
''')


@pytest.mark.slow
def test_report_identical_across_hash_seeds(tmp_path):
    from pathlib import Path

    fixture = Path(__file__).parent.parent / "fixtures" / "provenance_example.json"
    repo_root = Path(__file__).resolve().parents[2]
    script = tmp_path / "gen.py"
    script.write_text(_SCRIPT)
    outputs = []
    for seed in ("0", "777"):
        out = tmp_path / f"r_{seed}.html"
        env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(repo_root)}
        subprocess.run(
            [sys.executable, str(script), str(out), str(fixture)],
            check=True, env=env, cwd=repo_root, timeout=300,
        )
        text = out.read_text(encoding="utf-8")
        outputs.append(re.sub(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", "<TS>", text))
    assert outputs[0] == outputs[1]
    assert "Table truncated" in outputs[0]
