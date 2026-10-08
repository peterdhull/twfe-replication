"""Pairwise inference, design concentration screens, and marker semantics.

The stored paired-difference SE is computed directly from joint bootstrap draws.
It is retained rather than reconstructed by subtracting nearly equal variances.
The covariance identity is checked before any screen or figure is generated.
"""

import math

import numpy as np
import pandas as pd

KEY = ["family", "specification_id", "modern_estimator", "summary_period"]
SE_FLOOR = 1e-12


def boolean_series(values):
    """Accept actual booleans or their unambiguous CSV encodings."""
    if values.dtype == bool:
        return values.copy()
    mapped = values.map({True: True, False: False, "True": True, "False": False,
                         "true": True, "false": False, 1: True, 0: False})
    if mapped.isna().any():
        raise ValueError(f"Missing or invalid booleans in {values.name}")
    return mapped.astype(bool)


def concentration_failure(effective_clusters, maximum_share, rule, tolerance):
    """Strict thresholds with the documented floating-point boundary tolerance."""
    return ((effective_clusters < rule["effective_clusters_below"] - tolerance)
            | (maximum_share > rule["maximum_cluster_share_above"] + tolerance))


def paired_inference(frame, threshold=1.96):
    """Validate saved marginal/pair uncertainty and independently derive t ratios."""
    result = frame.copy()
    for column in ("twfe_estimate", "modern_estimate", "twfe_positive_se",
                   "modern_positive_se", "se_gap", "covariance_twfe_modern"):
        if not np.isfinite(result[column].to_numpy(float)).all():
            raise ValueError(f"Nonfinite values in {column}")
    sx = result.twfe_positive_se.to_numpy(float)
    sy = result.modern_positive_se.to_numpy(float)
    sg = result.se_gap.to_numpy(float)
    cov = result.covariance_twfe_modern.to_numpy(float)
    if (sx < 0).any() or (sy < 0).any() or (sg < 0).any():
        raise ValueError("Standard errors must be nonnegative")
    reconstructed = sx**2 + sy**2 - 2 * cov
    # Absolute tolerance covers CSV round-off and cancellation in identities.
    # Relative tolerance scales with the variances being subtracted, not the
    # potentially tiny difference variance.
    allowance = 1e-13 + 1e-8 * (sx**2 + sy**2 + 2 * np.abs(cov))
    if np.any(np.abs(sg**2 - reconstructed) > allowance):
        raise ValueError("Difference variance does not match marginal variances and covariance")
    if np.any(np.abs(cov) > sx * sy + allowance):
        raise ValueError("Invalid paired covariance")
    result["gap"] = result.twfe_estimate - result.modern_estimate
    result["modern_minus_twfe"] = -result.gap
    result["t_signed"] = np.divide(
        result.gap.to_numpy(float), sg, out=np.full(len(result), np.nan),
        where=sg > SE_FLOOR)
    result["abs_t"] = result.t_signed.abs()
    result["paired_normal_p_unadjusted"] = result.abs_t.map(
        lambda value: math.erfc(value / math.sqrt(2)) if np.isfinite(value) else np.nan)
    result["numerical_abs_t_above_1_96"] = result.abs_t > 1.96
    result["paired_difference_significant"] = result.abs_t > threshold
    return result


def screen_comparisons(frame, policy):
    """Recompute all inclusion/warning decisions from exported diagnostics."""
    if "summary_period" not in frame:
        frame = frame.assign(summary_period="static")
    if frame.duplicated(KEY).any():
        raise ValueError("Duplicate comparison keys")
    threshold = float(policy.get("pointwise_abs_t_threshold", 1.96))
    result = paired_inference(frame, threshold=threshold)
    for column in ("gap_structural_identity", "screen_prior_inference_warning"):
        result[column] = boolean_series(result[column])
    tolerance = float(policy["numerical_boundary_tolerance"])
    for role in ("twfe", "modern", "gap"):
        valid = ~result.gap_structural_identity if role == "gap" else np.ones(len(result), bool)
        geff = result.loc[valid, role + "_G_eff"]
        lmax = result.loc[valid, role + "_Lmax"]
        if not np.isfinite(geff).all() or not np.isfinite(lmax).all():
            raise ValueError(f"Nonfinite concentration metrics for {role}")
        if (geff < 1 - 1e-9).any() or not lmax.between(0, 1 + 1e-9).all():
            raise ValueError(f"Invalid concentration metrics for {role}")
    for role in ("twfe", "modern"):
        result[role + "_screen_exclude"] = concentration_failure(
            result[role + "_G_eff"], result[role + "_Lmax"],
            policy["hard_exclusion"], tolerance)
    result["screen_exclude"] = result.twfe_screen_exclude | result.modern_screen_exclude
    result["screen_leverage_warning"] = False
    for role in ("twfe", "modern", "gap"):
        flag = concentration_failure(result[role + "_G_eff"], result[role + "_Lmax"],
                                     policy["warning"], tolerance)
        if role == "gap":
            flag &= ~result.gap_structural_identity
        result[role + "_screen_warning"] = flag
        result["screen_leverage_warning"] |= flag
    result["screen_fragile"] = (result.screen_leverage_warning
                                | result.screen_prior_inference_warning)
    result["screen_status"] = np.where(
        result.screen_exclude, "excluded",
        np.where(result.screen_fragile, "retained_fragile", "retained"))
    result["screen_exclusion_reason"] = [
        "; ".join(role for role in ("twfe", "modern") if row[role + "_screen_exclude"])
        for row in result.to_dict("records")]
    result["fill_status"] = marker_status(result.abs_t, result.screen_fragile, threshold=threshold)
    result["sign_reversal"] = sign_reversal(result.twfe_estimate, result.modern_estimate)
    result["significant_reversal"] = result.sign_reversal & result.fill_status.eq("filled")
    return result


def sign_reversal(twfe, modern):
    """Strictly opposite nonzero signs, evaluated before coordinate truncation."""
    return ((twfe > 0) & (modern < 0)) | ((twfe < 0) & (modern > 0))


def marker_status(abs_t, fragile, threshold=1.96):
    """Filled iff pointwise paired |t| exceeds 1.96 and inference is not flagged."""
    return np.where(fragile, "fragile", np.where(abs_t > threshold, "filled", "hollow"))


def annotation_summary(frame):
    rows = []
    for (family, modern), group in frame.groupby(["family", "modern_estimator"]):
        reversal = sign_reversal(group.twfe_estimate, group.modern_estimate)
        filled = group.fill_status.eq("filled")
        rows.append(dict(family=family, modern_estimator=modern, comparisons=len(group),
                         pre_averages=int(group.summary_period.eq("pre").sum()),
                         post_averages=int(group.summary_period.eq("post").sum()),
                         scalar_summaries=int(group.summary_period.eq("static").sum()),
                         studies=int(group.study_id.nunique()),
                         fragile=int(group.screen_fragile.sum()),
                         sign_reversals=int(reversal.sum()),
                         sign_reversal_share=float(reversal.mean()),
                         filled=int(filled.sum()), filled_share=float(filled.mean()),
                         significant_reversals=int((reversal & filled).sum()),
                         significant_reversal_share=float((reversal & filled).mean())))
    return pd.DataFrame(rows)
