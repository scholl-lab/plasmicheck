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
}

# Heatmap colors matching Plotly config
HEATMAP_COLORS: dict[str, str] = {
    "not_contaminated": "lightblue",
    "unclear": "orange",
    "contaminated": "red",
}
