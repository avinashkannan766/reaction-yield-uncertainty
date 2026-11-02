"""Rebuild the benchmark table from the per-experiment summaries already saved under runs/.

    python -m chemyield.collect --dataset tdc_buchwald [--contains gkf] [--out results/tdc_buchwald_collected.md]

Use it when a long benchmark was interrupted (the benchmark only writes its table at the very end)
or to look at a subset of splits. Reads runs/<dataset>/<split>__<experiment>/summary.json.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _pm(v) -> str:
    return "nan" if v is None or v[0] != v[0] else f"{v[0]:.3f} ± {v[1]:.3f}"


def collect(dataset: str, root: str | Path = ".", contains: str | None = None) -> str:
    base = Path(root) / "runs" / dataset
    rows: dict[str, list[tuple]] = {}
    for f in sorted(base.glob("*__*/summary.json")):
        split, exp = f.parent.name.split("__", 1)
        if contains and contains not in split:
            continue
        s = json.loads(f.read_text())
        rows.setdefault(split, []).append((exp, s))
    lines = [f"# {dataset}: collected from runs/ (mean ± std over seeds / folds)", ""]
    for split, items in rows.items():
        lines += [f"## split: {split}", "",
                  "| experiment | R² | RMSE | 90% coverage | interval width |", "|---|---|---|---|---|"]
        for exp, s in items:
            lines.append(f"| {exp} | {_pm(s['r2'])} | {_pm(s['rmse'])} | {_pm(s['coverage'])} | "
                         f"{_pm(s['mean_width'])} |")
        lines.append("")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--contains", help="only splits whose name contains this (e.g. gkf)")
    ap.add_argument("--out", help="write markdown here (default: print)")
    args = ap.parse_args(argv)
    text = collect(args.dataset, contains=args.contains)
    if args.out:
        Path(args.out).write_text(text)
        print(f"Wrote {args.out}")
    else:
        print(text)


if __name__ == "__main__":
    main()
