"""View contrasts and optional descriptive mirror-tail screening."""

from fractions import Fraction

import numpy as np

from ._evaluate import _positive_int
from ._results import AuditResult, EvaluationResult


def _mirror_threshold(values, ratio):
    for threshold in np.unique(np.abs(values[values != 0])):
        positive = int(np.sum(values >= threshold))
        negative = int(np.sum(values <= -threshold))
        if positive and (1 + negative) / positive <= ratio:
            return float(threshold)
    return None


def _screening_summary(values, ratio, threshold):
    """Explain a decision without changing the threshold or evaluated scores."""
    positive_available = int(np.sum(values > 0))
    minimum = None
    if ratio is not None:
        fraction = Fraction(ratio)
        minimum = (fraction.denominator + fraction.numerator - 1) // fraction.numerator
        # Match floating-point comparisons used by the actual screening rule,
        # including the ordinary Python value 1 / 3.
        if minimum > 1 and 1 / (minimum - 1) <= ratio:
            minimum -= 1
    positive = int(np.sum(values >= threshold)) if threshold is not None else None
    negative = int(np.sum(values <= -threshold)) if threshold is not None else None
    if ratio is None:
        reason = "disabled"
        message = "Screening is disabled; continuous contrasts remain available."
    elif not len(values):
        reason = "no_eligible_groups"
        message = "No groups meet the minimum held-out support requirement."
    elif positive_available < minimum:
        reason = "insufficient_positive_groups"
        message = (
            f"Only {positive_available} eligible groups have positive contrast; "
            f"at least {minimum} are required at this mirror ratio, "
            "even with no negative-tail groups."
        )
    elif threshold is None:
        reason = "mirror_ratio_not_met"
        message = (
            "Enough positive groups are available, but no observed threshold "
            "meets the requested positive/negative tail ratio."
        )
    else:
        reason = "threshold_found"
        message = (
            f"Selected {positive} groups: (1 + {negative}) / {positive} "
            f"<= {ratio:g}. Candidates are descriptive, not confirmed artifacts."
        )
    return dict(
        reason=reason,
        message=message,
        eligible_count=len(values),
        positive_group_count=positive_available,
        minimum_positive_count=minimum,
        threshold=threshold,
        positive_tail_count=positive,
        negative_tail_count=negative,
        achieved_mirror_ratio=(1 + negative) / positive if positive else None,
        candidate_count=positive or 0,
    )


def audit(
    result,
    *,
    reference,
    contrast,
    groupby=None,
    min_test_cells=5,
    maximum_mirror_ratio=None,
):
    """Contrast group contributions from two views of one fixed grouping.

    Args:
        result (EvaluationResult): Output of evaluate; no model refitting occurs.
        reference (str): Baseline view, usually molecular.
        contrast (str): Compared view, usually QC. Difference is contrast-reference.
        groupby (str, optional): Required when multiple candidates were evaluated.
        min_test_cells (int): Minimum group count in EVERY test fold, default 5.
        maximum_mirror_ratio (float, optional): Enable descriptive screening with
            (1 + negative_tail_count) / positive_tail_count <= this value in (0,1].
            None (default) returns continuous contrasts without candidate calls.
            An exploratory starting point is 1 / 3; this is not an FDR level.

    Returns:
        AuditResult: One row per group, including ineligible groups with reasons.
            Differences are median-centered over eligible groups within this
            grouping ONLY. No pooling across independent samples is performed.
            metadata['screening'] explains the decision, tail counts and minimum
            positive-group requirement. Undefined threshold/tail values are None.

    Notes:
        Screening is not an FDR procedure or a test of significance. QC coupling
        does not establish a technical artifact or identify individual doublets.
    """
    if not isinstance(result, EvaluationResult):
        raise ValueError(
            "INVALID_RESULT: audit requires flat evaluation; inspect hierarchy conditional view scores separately."
        )
    _positive_int(min_test_cells, "min_test_cells")
    if maximum_mirror_ratio is not None:
        if (
            isinstance(maximum_mirror_ratio, bool)
            or not np.isfinite(maximum_mirror_ratio)
            or not 0 < maximum_mirror_ratio <= 1
        ):
            raise ValueError("INVALID_MIRROR_RATIO: expected a number in (0, 1].")
        maximum_mirror_ratio = float(maximum_mirror_ratio)
    available = result.table["groupby"].unique()
    if groupby is None:
        if len(available) != 1:
            raise ValueError(
                "AMBIGUOUS_GROUPING: specify groupby when auditing multiple candidates."
            )
        groupby = available[0]
    if groupby not in available:
        raise ValueError("MISSING_GROUPING: candidate was not evaluated.")
    if (
        reference == contrast
        or reference not in result.metadata["views"]
        or contrast not in result.metadata["views"]
    ):
        raise ValueError("INVALID_VIEW_PAIR: select two distinct evaluated views.")
    g = result.groups[result.groups["groupby"] == groupby]
    left = g[g.view == reference].set_index("group")
    right = g[g.view == contrast].set_index("group").reindex(left.index)
    frame = left[["n_cells", "proportion", "information_mass", "min_test_cells"]].copy()
    frame["reference_contribution"] = left.contribution
    frame["contrast_contribution"] = right.contribution
    frame["reference_recall"] = left.recall
    frame["contrast_recall"] = right.recall
    frame["difference"] = right.contribution - left.contribution
    frame["eligible"] = frame.min_test_cells >= min_test_cells
    frame["reason"] = np.where(
        frame.eligible,
        "",
        f"fewer than {min_test_cells} observations in at least one test fold",
    )
    center = (
        float(frame.loc[frame.eligible, "difference"].median())
        if frame.eligible.any()
        else None
    )
    frame["contrast"] = np.nan
    if center is not None:
        frame.loc[frame.eligible, "contrast"] = (
            frame.loc[frame.eligible, "difference"] - center
        )
    threshold = None
    if maximum_mirror_ratio is not None:
        threshold = _mirror_threshold(
            frame.loc[frame.eligible, "contrast"].to_numpy(), maximum_mirror_ratio
        )
        frame["candidate"] = (
            False
            if threshold is None
            else frame.eligible & (frame["contrast"] >= threshold)
        )
    assignments = result.predictions[
        (result.predictions["groupby"] == groupby)
        & (result.predictions.view == reference)
    ][["obs_name", "group"]].copy()
    metadata = dict(
        schema_version=1,
        groupby=groupby,
        reference=reference,
        contrast=contrast,
        min_test_cells=int(min_test_cells),
        maximum_mirror_ratio=maximum_mirror_ratio,
        center=center,
        threshold=threshold,
        centering_scope="eligible groups within one grouping",
        status="no_eligible_groups" if center is None else "ok",
        screening_status=(
            "disabled"
            if maximum_mirror_ratio is None
            else "no_threshold" if threshold is None else "threshold_found"
        ),
        screening=_screening_summary(
            frame.loc[frame.eligible, "contrast"].to_numpy(),
            maximum_mirror_ratio,
            threshold,
        ),
        interpretation="Group-level descriptive screening; no p-values or FDR guarantee",
        evaluation=result.metadata,
    )
    return AuditResult(frame.reset_index(), assignments, metadata)
