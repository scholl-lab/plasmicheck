"""Shared color constants matching Plotly default categorical colors."""

from __future__ import annotations

# Plotly's default categorical color sequence (plotly.colors.qualitative.Plotly)
PLOTLY_COLORS: list[str] = [
    "#636EFA",
    "#EF553B",
    "#00CC96",
    "#AB63FA",
    "#FFA15A",
    "#19D3F3",
    "#FF6692",
    "#B6E880",
    "#FF97FF",
    "#FECB52",
]

# Category-specific colors for assignment plots
ASSIGNMENT_COLORS: dict[str, str] = {
    "Plasmid": "#636EFA",
    "Human": "#EF553B",
    "Tied": "#00CC96",
    "Backbone_Only": "#D3D3D3",  # lightgray
    "Ambiguous": "#A9A9A9",  # darkgray
}

# Heatmap colors matching clinical status palette
HEATMAP_COLORS: dict[str, str] = {
    "not_contaminated": "#10b981",  # emerald-500
    "unclear": "#f59e0b",  # amber-500
    "contaminated": "#ef4444",  # red-500
}

# Verdict bar chart colors derived from heatmap colors
VERDICT_COLORS: dict[str, str] = {
    "Clean": HEATMAP_COLORS["not_contaminated"],
    "Unclear": HEATMAP_COLORS["unclear"],
    "Contaminated": HEATMAP_COLORS["contaminated"],
    "No Result": "#94a3b8",  # slate-400
}
