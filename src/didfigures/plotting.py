"""Standalone scatter renderer with optional public bibliographic labels."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

from .screening import sign_reversal

INK = "#172B3A"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK, "ytick.color": INK,
                     "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})


def render(data, studies, bounds, destination, scale, capped, source_pairs, policy,
           figure_config=None, named=False):
    """Write PNG, SVG, and PDF with independently computed screen/marker fields."""
    destination = Path(destination)
    figure_config = figure_config or {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    aligned = data.family.iloc[0] == "aligned_samples"
    study_meta = studies.set_index("study_id")
    papers = sorted(data.study_id.unique())
    if not set(papers).issubset(study_meta.index):
        raise ValueError("Missing study metadata")
    if named:
        required = ["short_title", "short_citation", "source_url"]
        if not set(required).issubset(study_meta.columns) or study_meta.loc[papers, required].isna().any().any():
            raise ValueError("Named figures require complete public study provenance")
        if not study_meta.loc[papers, "source_url"].str.startswith("https://doi.org/").all():
            raise ValueError("Named figure links must be official DOI URLs")
    colors = study_meta.color.to_dict()
    fig, axes = plt.subplots(1, 2, figsize=(14, 10.8) if named else (14, 9.6))
    fig.subplots_adjust(left=.066, right=.973, top=.875 if named else .895,
                        bottom=.375 if named else .305, wspace=.20)
    threshold = float(policy.get("pointwise_abs_t_threshold", 1.96))
    fig.text(.066, .975,
             f"{len(papers)} studies · Filled: paired |t| > {threshold:g} · Hollow: otherwise · Cross: fragile inference",
             fontsize=11, color="#4F5F6B")
    counts = data.drop_duplicates(["specification_id", "summary_period"]).groupby("study_id").size()
    plot_order = sorted(papers, key=lambda p: (-counts[p], study_meta.loc[p, "draw_order_tiebreak"]))
    lo, hi = bounds
    cap_text = f"±{hi:g}" if np.isclose(lo, -hi) else f"[{lo:g}, {hi:g}]"
    axis_units = "TWFE SEs" if scale == "twfe-se" else "outcome SDs"
    checks, annotations = [], []
    for ax, modern in zip(axes, ("csa", "bjs")):
        rows = data.loc[data.modern_estimator.eq(modern)]
        visible = rows.plot_x.between(lo, hi) & rows.plot_y.between(lo, hi)
        shown = rows.copy()
        shown["capped"] = ~visible if capped else False
        shown["display_x"] = shown.plot_x.clip(lo, hi) if capped else shown.plot_x
        shown["display_y"] = shown.plot_y.clip(lo, hi) if capped else shown.plot_y
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        ax.set_axisbelow(True)
        ax.grid(color="#E5E9ED", lw=.6)
        ax.axhline(0, color="#BEC6CD", lw=.85, zorder=1)
        ax.axvline(0, color="#BEC6CD", lw=.85, zorder=1)
        ax.plot([lo, hi], [lo, hi], color="#343E47", lw=1.1,
                linestyle=(0, (5, 4)), zorder=2)
        point_count = 0
        for study in plot_order:
            group = shown.loc[shown.study_id.eq(study)]
            for (truncated, status, period), subgroup in group.groupby(["capped", "fill_status", "summary_period"]):
                color = colors[study]
                is_event = period in ("pre", "post")
                if is_event and (status == "fragile" or truncated):
                    # An outline preserves the event-summary shape behind a cross;
                    # an additional diamond identifies clipping without hiding it.
                    if truncated:
                        ax.scatter(subgroup.display_x, subgroup.display_y, s=91,
                                   facecolors="none", edgecolors=color, marker="D",
                                   linewidths=.8, clip_on=False, zorder=6)
                    if status == "fragile":
                        ax.scatter(subgroup.display_x, subgroup.display_y, s=45,
                                   facecolors="none", edgecolors=color, marker="s",
                                   linewidths=.8, clip_on=not truncated, zorder=6)
                if status == "fragile":
                    collection = ax.scatter(subgroup.display_x, subgroup.display_y,
                                            s=39, c=color, marker="x", linewidths=1.25,
                                            clip_on=not truncated, zorder=7)
                    if truncated and not is_event:
                        ax.scatter(subgroup.display_x, subgroup.display_y, s=55,
                                   facecolors="none", edgecolors=color, marker="D",
                                   linewidths=.9, clip_on=False, zorder=6)
                else:
                    filled = status == "filled"
                    collection = ax.scatter(subgroup.display_x, subgroup.display_y,
                                            s=37 if truncated or is_event else 31,
                                            facecolors=color if filled else "none",
                                            edgecolors=color, marker="s" if is_event else ("D" if truncated else "o"),
                                            linewidths=.8 if filled else 1.0, alpha=.90,
                                            clip_on=not truncated, zorder=(8 if is_event else (5 if filled else 3)))
                    assert (len(collection.get_facecolors()) > 0) == filled
                np.testing.assert_array_equal(collection.get_offsets(),
                                              subgroup[["display_x", "display_y"]].to_numpy())
                point_count += len(collection.get_offsets())
        assert point_count == len(shown)
        ax.set_title(f"{modern.upper()} vs TWFE", loc="left", fontsize=14, fontweight="bold", pad=28)
        detail = f"{len(rows)} retained / {source_pairs[modern]} comparisons"
        if capped:
            detail += f"; {int((~visible).sum())} capped (diamonds)"
        ax.text(0, 1.025, detail, transform=ax.transAxes, fontsize=9.5, color="#5A6873")
        reversals = sign_reversal(rows.twfe_estimate, rows.modern_estimate)
        reversal_n = int(reversals.sum())
        filled_n = int(shown.fill_status.eq("filled").sum())
        significant_reversal_n = int((reversals & shown.fill_status.eq("filled")).sum())
        label = (f"Sign-reversals: {100*reversal_n/len(rows):.1f}%\n"
                 f"Significant differences: {100*filled_n/len(rows):.1f}%\n"
                 f"Significant + Reversal: {100*significant_reversal_n/len(rows):.1f}%")
        annotation = ax.text(.035 if capped else .965, .975 if capped else .045,
                             label, transform=ax.transAxes,
                             ha="left" if capped else "right",
                             va="top" if capped else "bottom", fontsize=10.2,
                             linespacing=1.55, color=INK, zorder=9,
                             bbox=dict(boxstyle="round,pad=.55", facecolor="white",
                                       edgecolor="#CBD3D9", linewidth=.7, alpha=.97))
        annotations.append((annotation, ax, shown[["display_x", "display_y"]].to_numpy()))
        ax.set_xlabel(("Sample-aligned TWFE" if aligned else "TWFE") + f" estimate ({axis_units})",
                      fontsize=11, labelpad=9)
        ax.set_ylabel(f"{modern.upper()} estimate ({axis_units})", fontsize=11, labelpad=8)
        if capped:
            ax.set_xticks(np.linspace(lo, hi, 5))
            ax.set_yticks(np.linspace(lo, hi, 5))
        else:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=7))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
        for spine in ax.spines.values():
            spine.set_color("#A0ABB4")
        checks.append(dict(estimator=modern, total=len(rows), plotted=point_count,
                           capped=int(shown.capped.sum()), uncapped=int((~shown.capped).sum()),
                           bounds=list(bounds), fill_counts=shown.fill_status.value_counts().to_dict(),
                           sign_reversal_count=reversal_n, sign_reversal_share=reversal_n/len(rows),
                           filled_count=filled_n, filled_share=filled_n/len(rows),
                           significant_reversal_count=significant_reversal_n,
                           significant_reversal_share=significant_reversal_n/len(rows),
                           study_count=int(rows.study_id.nunique()),
                           pre_average_count=int(rows.summary_period.eq("pre").sum()),
                           post_average_count=int(rows.summary_period.eq("post").sum()),
                           scalar_count=int(rows.summary_period.eq("static").sum()),
                           visible_summary_points=list(zip(shown.specification_id, shown.summary_period)),
                           visible_specifications=shown.specification_id.tolist()))

    order = sorted(papers, key=lambda p: (study_meta.loc[p, "publication_year"], p))
    legend_artists = []
    columns = int(figure_config.get("legend_columns", 4))
    if columns < 1:
        raise ValueError("legend_columns must be positive")
    quotient, remainder = divmod(len(order), columns)
    column_counts = [quotient + int(col < remainder) for col in range(columns)]
    legend_rows = max(column_counts)
    positions = [(col, row) for col, count in enumerate(column_counts) for row in range(count)]
    legend_top = (.291 if legend_rows <= 6 else .307) if named else (.235 if legend_rows <= 6 else .251)
    legend_step = (legend_top - .143) / max(1, legend_rows - 1)
    for study, (col, row) in zip(order, positions):
        x, y = .077 + col * (.924 / columns), legend_top - row * legend_step
        dot = Line2D([x], [y+.001], marker="o", linestyle="", markersize=5.5,
                     markerfacecolor=colors[study], markeredgecolor="none", transform=fig.transFigure)
        fig.add_artist(dot)
        if named:
            title = fig.text(x+.012, y+.003, study_meta.loc[study, "short_title"],
                             fontsize=8.4, va="bottom")
            citation = fig.text(x+.012, y-.009, study_meta.loc[study, "short_citation"],
                                fontsize=7.2, color="#61717D", va="bottom")
            citation.set_url(study_meta.loc[study, "source_url"])
            legend_artists.extend([title, citation])
        else:
            legend_artists.append(fig.text(x+.012, y-.002,
                f"Study {study} ({int(study_meta.loc[study, 'publication_year'])})", fontsize=9.5, va="center"))
    hard, warning = policy["hard_exclusion"], policy["warning"]
    notes = [
        "Notes: Circles: static specifications. Squares: event-study averages, with pre and post entering separately. Pre averages are diagnostics, not treatment effects.",
        "Shares use retained comparisons, including flagged points. Significant + Reversal: opposite nonzero signs and a filled marker.",
        (f"Diamonds: coordinates capped at {cap_text}; displayed gaps may differ from actual gaps." if capped else
         "Event studies contribute separate pre- and post-period averages, never individual horizons."),
        (f"Excluded if either estimate has effective clusters <{hard['effective_clusters_below']:g} "
         f"or maximum cluster leverage >{100*hard['maximum_cluster_share_above']:g}%."),
        (f"Crosses: effective clusters <{warning['effective_clusters_below']:g} "
         f"or maximum leverage >{100*warning['maximum_cluster_share_above']:g}% for either estimate "
         "or their difference, or prior inference warnings; significance fill suppressed."),
        "Leverage uses squared outcome-linear weights by the existing bootstrap clusters; it is a design diagnostic, not a variance guarantee.",
        ("Exact-union sample alignment does not equate estimator weights or estimands." if aligned else
         "TWFE retains its eligible sample; modern estimators use their own supported targets."),
        ("Both axes divide by this pair's saved TWFE SE; the plotted gap is not the paired-difference t statistic. Common horizons are unchanged. Fill uses paired-bootstrap normal tests, without multiplicity adjustment. Dashed line: equality."
         if scale == "twfe-se" else
         "Original-SD scaling and common-horizon event averages are unchanged. Fill uses existing paired-bootstrap normal tests, without multiplicity adjustment. Dashed line: equality."),
    ]
    if 8 in papers:
        notes.append("Study 8 uses unrestricted TWFE; its restricted published fit is retained separately.")
    if 4 in papers:
        notes.append("Study 4 clusters 143 builders despite only three treatment-rollout regions.")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    footer_font = FontProperties(family=plt.rcParams["font.family"], size=8.3)
    footer_left, footer_right = .066, .973
    footer_width = (footer_right - footer_left) * fig.bbox.width

    def text_width(value):
        return renderer.get_text_width_height_descent(value, footer_font, ismath=False)[0]

    lines, current = [], []
    for word in " ".join(notes).split():
        if current and text_width(" ".join(current + [word])) > footer_width:
            lines.append(current)
            current = []
        current.append(word)
    if current:
        lines.append(current)
    footer_artists = []
    footer_top = .117
    footer_step = min(.017, (.117 - .02) / max(1, len(lines) - 1))
    assert footer_step >= .012, "Figure notes require too many lines"
    for i, words in enumerate(lines):
        y = footer_top - footer_step*i
        assert y >= .019 - 1e-12
        if i == len(lines)-1 or len(words) == 1:
            footer_artists.append(fig.text(footer_left, y, " ".join(words),
                                           fontproperties=footer_font, color="#4F5F6B"))
            continue
        widths = [text_width(word) for word in words]
        gap = (footer_width - sum(widths)) / (len(words)-1)
        assert gap >= text_width(" ") - .1
        x = footer_left * fig.bbox.width
        for word, width in zip(words, widths):
            footer_artists.append(fig.text(x / fig.bbox.width, y, word,
                                           fontproperties=footer_font, color="#4F5F6B"))
            x += width + gap
    fig.canvas.draw()
    for artist in footer_artists:
        bbox = artist.get_window_extent(renderer)
        assert bbox.x0 >= footer_left*fig.bbox.width - 1
        assert bbox.x1 <= footer_right*fig.bbox.width + 1
        assert bbox.y0 >= 0
        assert not any(bbox.overlaps(a.get_window_extent(renderer)) for a in legend_artists)
    for artist in legend_artists:
        bbox = artist.get_window_extent(renderer)
        assert 0 <= bbox.x0 <= bbox.x1 <= fig.bbox.width
        assert 0 <= bbox.y0 <= bbox.y1 <= fig.bbox.height
        assert not any(bbox.overlaps(ax.xaxis.label.get_window_extent(renderer)) for ax in axes)
    if named:
        boxes = [artist.get_window_extent(renderer) for artist in legend_artists]
        assert not any(box.overlaps(other) for i, box in enumerate(boxes) for other in boxes[i+1:]), "Named legend labels overlap"
    for ax in axes:
        transform = ax.transData.transform
        np.testing.assert_allclose(np.linalg.norm(transform([1, 0])-transform([0, 0])),
                                   np.linalg.norm(transform([0, 1])-transform([0, 0])), rtol=1e-12)
    for annotation, ax, xy in annotations:
        screens = ax.transData.transform(xy)
        clearance = 5 * fig.dpi / 72
        default_position = (.035, .79) if scale == "twfe-se" else None
        position = figure_config.get("scales", {}).get(scale, {}).get("capped_annotation_axes", default_position)
        if capped and position is not None:
            candidates = [(tuple(position), "left", "top")]
        else:
            candidates = [(annotation.get_position(), annotation.get_ha(), annotation.get_va())]
            candidates += [((x, y), "left", "bottom")
                           for y in np.linspace(.80, .035, 28) for x in np.linspace(.035, .60, 24)]
        placed = False
        for position, ha, va in candidates:
            annotation.set_position(position)
            annotation.set_ha(ha)
            annotation.set_va(va)
            annotation.update_bbox_position_size(renderer)
            bbox = annotation.get_bbox_patch().get_window_extent(renderer)
            if not (ax.bbox.contains(bbox.x0, bbox.y0) and ax.bbox.contains(bbox.x1, bbox.y1)):
                continue
            if any(bbox.padded(clearance).contains(x, y) for x, y in screens):
                continue
            placed = True
            break
        assert placed, "No point-free annotation position found"
        checks[axes.tolist().index(ax)]["annotation_position_axes"] = list(annotation.get_position())
        checks[axes.tolist().index(ax)]["annotation_point_clearance_pt"] = 5
    fig.canvas.draw()
    for extension in figure_config.get("formats", ("png", "svg", "pdf")):
        if extension not in ("png", "svg", "pdf"):
            raise ValueError(f"Unsupported output format: {extension}")
        fig.savefig(destination.with_suffix("." + extension), dpi=190, facecolor="white")
    links = {a.get_url() for a in fig.findobj() if hasattr(a, "get_url") and a.get_url()}
    if named:
        assert links.issubset(set(study_meta.loc[papers, "source_url"]))
    else:
        assert not links
    plt.close(fig)
    return checks
