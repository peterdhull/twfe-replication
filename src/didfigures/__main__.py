"""Command-line entry point for self-contained figure replication."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from . import __version__
from .plotting import render
from .screening import annotation_summary, screen_comparisons


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_clean(value):
    if isinstance(value, dict):
        return {str(k): json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(v) for v in value]
    if isinstance(value, np.generic):
        return json_clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def run(root, output, scale="both", view="both", summaries="all"):
    root = Path(root).resolve()
    output = Path(output)
    if not output.is_absolute():
        output = root / output
    output.mkdir(parents=True, exist_ok=True)
    input_files = [root / "data/results/comparisons.csv", root / "data/results/studies.csv",
                   root / "config/screen.json", root / "config/figures.json",
                   root / "data/results/pre_comparisons.csv"]
    protected = {str(p.relative_to(root)).replace("\\", "/"): sha256(p) for p in input_files}
    data = pd.read_csv(input_files[0], float_precision="round_trip")
    pre = pd.read_csv(input_files[4], float_precision="round_trip")
    if not pre.summary_period.eq("pre").all():
        raise ValueError("pre_comparisons.csv must contain only pooled pre summaries")
    data = pd.concat([data, pre], ignore_index=True)
    if summaries == "treatment":
        data = data.loc[~data.summary_period.eq("pre")].copy()
    elif summaries == "pre":
        data = data.loc[data.summary_period.eq("pre")].copy()
    elif summaries != "all":
        raise ValueError("Unknown summary selection")
    studies = pd.read_csv(input_files[1])
    policy = json.loads(input_files[2].read_text(encoding="utf-8-sig"))
    figure_config = json.loads(input_files[3].read_text(encoding="utf-8-sig"))
    if set(data.family) != {"baseline", "aligned_samples"}:
        raise ValueError("Expected baseline and aligned_samples families")
    if studies.study_id.duplicated().any():
        raise ValueError("Duplicate study metadata")
    decisions = screen_comparisons(data, policy)
    retained = decisions.loc[~decisions.screen_exclude].copy()
    summaries = annotation_summary(retained)
    period_summaries = pd.concat([
        annotation_summary(group).assign(summary_period=period)
        for period, group in retained.groupby("summary_period")], ignore_index=True)
    decisions.to_csv(output / "screening_decisions.csv", index=False, float_format="%.17g")
    retained.to_csv(output / "retained_plot_data.csv", index=False, float_format="%.17g")
    summaries.to_csv(output / "annotation_summary.csv", index=False, float_format="%.17g")
    period_summaries.to_csv(output / "annotation_summary_by_period.csv", index=False, float_format="%.17g")
    scales = ("twfe-se", "outcome-sd") if scale == "both" else (scale,)
    views = ("capped", "full") if view == "both" else (view,)
    figures, ranges = {}, {}
    coordinates = []
    for unit in scales:
        displayed = retained.copy()
        if unit == "twfe-se":
            divisor = displayed.twfe_positive_se.to_numpy(float)
            if not np.isfinite(divisor).all() or (divisor <= 0).any():
                raise ValueError("TWFE-SE presentation requires strictly positive finite TWFE SEs")
        else:
            divisor = np.ones(len(displayed))
        capped_bounds = tuple(float(v) for v in figure_config["scales"][unit]["capped_bounds"])
        if len(capped_bounds) != 2 or not capped_bounds[0] < capped_bounds[1]:
            raise ValueError("Expected two increasing capped_bounds")
        minpad = max(abs(v) for v in capped_bounds) / 10
        displayed["plot_x"] = displayed.twfe_estimate / divisor
        displayed["plot_y"] = displayed.modern_estimate / divisor
        displayed["axis_scale"] = unit
        coordinates.append(displayed[["family", "specification_id", "modern_estimator", "summary_period", "axis_scale",
                                      "plot_x", "plot_y"]])
        values = displayed[["plot_x", "plot_y"]].to_numpy().ravel()
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite plotting coordinates")
        pad = max(minpad, .04 * (values.max() - values.min()))
        full_bounds = (min(capped_bounds[0], float(values.min())-pad),
                       max(capped_bounds[1], float(values.max())+pad))
        ranges[unit] = dict(capped=list(capped_bounds), full=list(full_bounds))
        figures[unit] = {}
        for family in ("baseline", "aligned_samples"):
            frame = displayed.loc[displayed.family.eq(family)]
            source_pairs = data.loc[data.family.eq(family), "modern_estimator"].value_counts().to_dict()
            figures[unit][family] = {}
            for perspective in views:
                is_capped = perspective == "capped"
                destination = output / unit / f"{family}_{perspective}"
                checks = render(frame, studies, capped_bounds if is_capped else full_bounds,
                                destination, unit, is_capped, source_pairs, policy, figure_config)
                figures[unit][family][perspective] = checks
                for check in checks:
                    expected = summaries.loc[(summaries.family == family)
                                             & (summaries.modern_estimator == check["estimator"])].iloc[0]
                    assert check["total"] == expected.comparisons
                    assert check["filled_count"] == expected.filled
                    assert check["sign_reversal_count"] == expected.sign_reversals
                    assert check["significant_reversal_count"] == expected.significant_reversals
    pd.concat(coordinates, ignore_index=True).to_csv(output / "plot_coordinates.csv", index=False,
                                                   float_format="%.17g")
    for file, digest in protected.items():
        assert sha256(root / file) == digest, f"Input changed: {file}"
    counts = []
    for (family, modern), group in decisions.groupby(["family", "modern_estimator"]):
        counts.append(dict(family=family, modern_estimator=modern, source=len(group),
                           retained=int((~group.screen_exclude).sum()),
                           excluded=int(group.screen_exclude.sum())))
    report = dict(status="PASS", package_version=__version__,
                  software={"numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
                  input_sha256=protected, source_studies=int(data.study_id.nunique()),
                  source_specifications=int(data.specification_id.nunique()),
                  source_summary_points=len(data.drop_duplicates(["specification_id", "summary_period"])),
                  retained_studies=int(retained.study_id.nunique()),
                  retained_specifications=int(retained.specification_id.nunique()),
                  retained_summary_points=len(retained.drop_duplicates(["specification_id", "summary_period"])),
                  source_pairs=len(data), source_and_retained_counts=counts,
                  annotation_summary=summaries.to_dict("records"),
                  annotation_summary_by_period=period_summaries.to_dict("records"), bounds=ranges, figures=figures,
                  source_screen_and_significance_flags_recomputed=True,
                  paired_variance_identity_checked=True,
                  direct_difference_se_used_to_avoid_covariance_cancellation=True,
                  inputs_unchanged=True)
    (output / "validation.json").write_text(json.dumps(json_clean(report), indent=2, allow_nan=False),
                                            encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(output),
                      "source_and_retained_counts": counts}, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2],
                        help="Repository root (default: this source checkout)")
    parser.add_argument("--output", type=Path, default=Path("outputs"), help="Output directory, relative to root unless absolute")
    parser.add_argument("--scale", choices=("twfe-se", "outcome-sd", "both"), default="both")
    parser.add_argument("--view", choices=("capped", "full", "both"), default="both")
    parser.add_argument("--summaries", choices=("all", "treatment", "pre"), default="all",
                        help="Include scalar/post and separate pre averages (default), or a subset")
    arguments = parser.parse_args()
    run(arguments.root, arguments.output, arguments.scale, arguments.view, arguments.summaries)


if __name__ == "__main__":
    main()
