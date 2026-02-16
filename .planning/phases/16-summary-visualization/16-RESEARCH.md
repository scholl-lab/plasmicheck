# Phase 16: Summary Visualization - Research

**Researched:** 2026-02-16
**Domain:** Multi-sample summary report visualization (Plotly + Matplotlib)
**Confidence:** HIGH

## Summary

This phase adds two new stacked bar charts to the multi-sample summary reports: one grouped by plasmid showing verdict counts, and one grouped by sample. The implementation follows the exact pattern established in Phase 13 (extended metrics heatmaps) — Plotly interactive charts with Matplotlib PNG fallbacks, integrated into the existing `generate_summary_reports.py` and `summary_template.html`.

Plotly stacked bar charts use `go.Bar()` traces with `barmode='stack'` to layer verdict categories (Contaminated, Unclear, Clean, No Result) as colored segments on each bar. Matplotlib uses the `bottom` parameter to stack bars iteratively. Both support color customization, hover text, and x-axis label rotation for readability.

**Primary recommendation:** Use Plotly `graph_objects.Bar()` with `customdata` for hover lists, aggregate verdict data with pandas `groupby().size().unstack().fillna(0)`, and create Matplotlib equivalents using iterative `ax.bar()` with cumulative `bottom` arrays.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| plotly | ~5.x | Interactive stacked bar charts | Already used for all summary charts (heatmaps, boxplots) |
| matplotlib | ~3.x | Static PNG fallbacks | Already used for `--static-report` mode |
| pandas | ~2.x | Data aggregation and reshaping | Already used throughout PlasmiCheck |
| jinja2 | ~3.x | HTML template rendering | Already used for report generation |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | ~1.x | Array operations for matplotlib stacking | Calculating cumulative `bottom` values |
| seaborn | ~0.x | Optional matplotlib styling | Already imported in matplotlib_backend.py |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| graph_objects.Bar() | plotly.express.bar() | px.bar() simpler but less control over customdata/hover; go.Bar() needed for complex hover lists |
| Stacked bars | Grouped bars (barmode='group') | Grouped shows absolute counts side-by-side but doesn't show total; stacked better for composition view |

**Installation:**
No new dependencies — all libraries already in `pyproject.toml` dev requirements.

## Architecture Patterns

### Recommended Integration Structure
```
plasmicheck/scripts/
├── generate_summary_reports.py    # Add verdict bar chart functions here
│   ├── plot_verdict_bar_by_plasmid()   # New function
│   ├── plot_verdict_bar_by_sample()    # New function
│   └── generate_report()               # Update to include new charts
├── plotting/
│   ├── matplotlib_backend.py      # Add matplotlib fallback functions
│   │   ├── generate_verdict_bar_matplotlib()
│   └── colors.py                  # Already has verdict colors
└── templates/
    └── summary_template.html      # Add new section at top
```

### Pattern 1: Verdict Data Aggregation (Pandas)
**What:** Transform summary_df (Sample, Plasmid, Category, Value) into counts by verdict
**When to use:** Before plotting, to prepare data for stacked bars
**Example:**
```python
# Source: pandas stacked bar chart documentation + existing PlasmiCheck pattern
import pandas as pd

# Filter to verdict rows
verdict_df = summary_df[summary_df["Category"] == "Verdict"].copy()

# Map full verdict strings to short labels
verdict_map = {
    "Sample is contaminated with plasmid DNA": "Contaminated",
    "Sample contamination status is unclear": "Unclear",
    "Sample is not contaminated with plasmid DNA": "Clean",
}
verdict_df["VerdictShort"] = verdict_df["Value"].map(verdict_map)

# Aggregate by plasmid: count verdicts per plasmid
by_plasmid = (
    verdict_df.groupby(["Plasmid", "VerdictShort"])
    .size()
    .unstack(fill_value=0)  # Rows=plasmids, cols=verdicts, fill missing with 0
)

# Aggregate by sample: count verdicts per sample
by_sample = (
    verdict_df.groupby(["Sample", "VerdictShort"])
    .size()
    .unstack(fill_value=0)
)

# Handle missing verdicts (sample-plasmid pairs with no result)
# Add "No Result" category by checking all_pairs - existing_pairs
# This requires cross-product of samples x plasmids, then left join
```

### Pattern 2: Plotly Stacked Bar with Custom Hover
**What:** Create interactive stacked bar chart with hover showing verdict + count + list of names
**When to use:** For interactive HTML reports with detailed hover info
**Example:**
```python
# Source: Plotly Bar charts documentation + hover text formatting
import plotly.graph_objects as go

# Assume aggregated_df has columns: [Contaminated, Unclear, Clean, No Result]
# Index is plasmid names or sample names

fig = go.Figure()

# Define verdict categories and colors (matching existing config)
verdict_categories = ["Clean", "Unclear", "Contaminated", "No Result"]
verdict_colors = {
    "Clean": "lightblue",      # From config.json plot_summary.colors.not_contaminated
    "Unclear": "orange",        # From config.json plot_summary.colors.unclear
    "Contaminated": "red",      # From config.json plot_summary.colors.contaminated
    "No Result": "#D3D3D3",     # Gray for missing data
}

# Build list of sample/plasmid names for each verdict (for hover)
# Example: hover_data["Contaminated"]["plasmid1"] = ["sample1", "sample2", ...]
hover_data = {}  # verdict -> bar_name -> list of names

# Add one trace per verdict category (stacks automatically)
for verdict in verdict_categories:
    counts = aggregated_df[verdict]  # Series with bar names as index

    # Build hover text with name lists (cap at 10, show "...and N more")
    hover_texts = []
    for bar_name in aggregated_df.index:
        name_list = hover_data[verdict].get(bar_name, [])
        if len(name_list) <= 10:
            names_str = "<br>".join(name_list)
        else:
            names_str = "<br>".join(name_list[:10]) + f"<br>...and {len(name_list)-10} more"
        hover_texts.append(
            f"<b>{verdict}</b><br>Count: {counts[bar_name]}<br>{names_str}"
        )

    fig.add_trace(go.Bar(
        x=aggregated_df.index,
        y=counts,
        name=verdict,
        marker_color=verdict_colors[verdict],
        hovertemplate="%{customdata}<extra></extra>",
        customdata=hover_texts,
    ))

fig.update_layout(
    barmode='stack',
    xaxis_title="Plasmid",
    yaxis_title="Count",
    title="Verdicts by Plasmid",
    xaxis=dict(tickangle=45),  # Rotate labels for readability
    height=600,
)

fig.write_html(output_path, include_plotlyjs=False, full_html=False)
```

### Pattern 3: Matplotlib Stacked Bar Fallback
**What:** Generate static PNG with stacked bars using iterative bottom parameter
**When to use:** When `--static-report` flag is set
**Example:**
```python
# Source: Matplotlib stacked bar chart documentation
import matplotlib.pyplot as plt
import numpy as np

# Assume aggregated_df has columns: [Clean, Unclear, Contaminated, No Result]
# Index is plasmid names or sample names

verdict_categories = ["Clean", "Unclear", "Contaminated", "No Result"]
verdict_colors = {
    "Clean": "lightblue",
    "Unclear": "orange",
    "Contaminated": "red",
    "No Result": "#D3D3D3",
}

fig, ax = plt.subplots(figsize=(10, 6))

# Initialize bottom at zero for first layer
bottom = np.zeros(len(aggregated_df))

# Stack each verdict category
for verdict in verdict_categories:
    counts = aggregated_df[verdict].values
    ax.bar(
        aggregated_df.index,
        counts,
        bottom=bottom,
        label=verdict,
        color=verdict_colors[verdict],
        width=0.6,
    )
    bottom += counts  # Update bottom for next layer

ax.set_xlabel("Plasmid")
ax.set_ylabel("Count")
ax.set_title("Verdicts by Plasmid")
ax.legend()

# Rotate x-axis labels if many items
if len(aggregated_df) > 5:
    plt.xticks(rotation=45, ha='right')

plt.tight_layout()
fig.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Pattern 4: Template Integration (Side-by-Side Charts)
**What:** Display two charts side by side at top of summary report
**When to use:** Integrating new charts into existing HTML template
**Example:**
```html
<!-- Source: Bootstrap grid pattern (already used in summary_template.html) -->
<!-- Add this section BEFORE "Contamination Analysis" card -->

<div class="card mb-3">
    <div class="card-header" id="headingVerdictOverview">
        <h5 class="mb-0">
            <button class="btn btn-link" data-toggle="collapse" data-target="#collapseVerdictOverview" aria-expanded="true">
                Verdict Overview
            </button>
        </h5>
    </div>
    <div id="collapseVerdictOverview" class="collapse show" aria-labelledby="headingVerdictOverview">
        <div class="card-body">
            <div class="row">
                <div class="col-md-6">
                    <h3>Verdicts by Plasmid</h3>
                    <div class="plot">
                        {% if interactive %}
                            {{ verdict_by_plasmid_content|safe }}
                        {% else %}
                            <img src="{{ verdict_by_plasmid_png }}" alt="Verdicts by Plasmid">
                        {% endif %}
                    </div>
                </div>
                <div class="col-md-6">
                    <h3>Verdicts by Sample</h3>
                    <div class="plot">
                        {% if interactive %}
                            {{ verdict_by_sample_content|safe }}
                        {% else %}
                            <img src="{{ verdict_by_sample_png }}" alt="Verdicts by Sample">
                        {% endif %}
                    </div>
                </div>
            </div>
        </div>
    </div>
</div>
```

### Anti-Patterns to Avoid
- **Using plotly.express for complex hover:** px.bar() doesn't support customdata lists well; use go.Bar() instead
- **Forgetting fillna(0):** Missing verdicts will create NaN values; always use `.unstack(fill_value=0)`
- **Hardcoding colors:** Use existing `HEATMAP_COLORS` from `colors.py` for consistency
- **Ignoring "No Result":** Must handle sample-plasmid pairs with no verdict (e.g., failed runs)

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Building hover text with lists | Manual string concatenation in loop | Plotly `customdata` + `hovertemplate` | Official pattern for variable-length hover content |
| Rotating x-axis labels | CSS transforms or manual plotting | Plotly `xaxis.tickangle` / Matplotlib `plt.xticks(rotation=45)` | Built-in, handles edge cases (label collision) |
| Stacking bars in matplotlib | Manual y-coordinate calculation | `bottom` parameter + cumulative sum | Matplotlib's designed approach, handles negative values |
| Missing data handling | Custom logic to detect missing pairs | Pandas `unstack(fill_value=0)` | Automatically aligns indices and fills gaps |
| Color consistency | Separate color definitions | Import from `colors.py` | Single source of truth already exists |

**Key insight:** Plotly's `customdata` and `hovertemplate` are designed for exactly this use case (dynamic lists in hover), and pandas `unstack()` automatically handles the cross-product logic for missing data.

## Common Pitfalls

### Pitfall 1: Hover Lists Exceeding Browser Limits
**What goes wrong:** With 50+ samples/plasmids, hover tooltips become unreadable
**Why it happens:** No cap on list length → massive hover boxes
**How to avoid:** Cap list display at ~10 items, show "...and N more" for remainder
**Warning signs:** Hover tooltips extend off-screen, slow hover rendering

### Pitfall 2: Inconsistent Verdict Color Mapping
**What goes wrong:** Charts use different colors than heatmaps/tables
**Why it happens:** Hardcoded colors instead of using config.json
**How to avoid:** Import `HEATMAP_COLORS` from `plotting/colors.py`
**Warning signs:** User confusion about verdict meanings, inconsistent legend

### Pitfall 3: X-Axis Label Overlap with Many Bars
**What goes wrong:** Bar labels overlap and become unreadable
**Why it happens:** Fixed chart width with many categories
**How to avoid:** Set `xaxis.tickangle=45` for >5 items; for 30+, increase chart width proportionally
**Warning signs:** Labels rendering on top of each other, user complaints

### Pitfall 4: "No Result" Not Distinguished from "Clean"
**What goes wrong:** Users can't tell missing data from actual clean verdicts
**Why it happens:** Missing verdicts filled with 0 but not shown as separate category
**How to avoid:** Add "No Result" as 4th verdict category with distinct gray color
**Warning signs:** Confusion about sample-plasmid pairs with no report

### Pitfall 5: matplotlib bottom Array Mutation
**What goes wrong:** Stacked bars render incorrectly (heights wrong)
**Why it happens:** Reusing same `bottom` array without copy, or in-place modification
**How to avoid:** Use `bottom += counts` (creates new array) or `bottom = bottom + counts`
**Warning signs:** Second/third bars starting at wrong y-position

### Pitfall 6: Plotly JS Not Loaded for New Charts
**What goes wrong:** Interactive charts fail to render in HTML
**Why it happens:** Forgetting to pass plotly_mode/version to template
**How to avoid:** Follow exact pattern from existing charts (heatmap, boxplot) — they already work
**Warning signs:** Empty divs in HTML, browser console errors about Plotly undefined

## Code Examples

Verified patterns from official sources:

### Aggregating Verdict Data with Missing Pairs
```python
# Source: pandas groupby + PlasmiCheck pattern for summary_df
import pandas as pd

# Load verdict data from summary_df
verdict_df = summary_df[summary_df["Category"] == "Verdict"].copy()

# Map full verdict strings to short labels
verdict_map = {
    "Sample is contaminated with plasmid DNA": "Contaminated",
    "Sample contamination status is unclear": "Unclear",
    "Sample is not contaminated with plasmid DNA": "Clean",
}
verdict_df["VerdictShort"] = verdict_df["Value"].map(verdict_map)

# Get all unique samples and plasmids
all_samples = summary_df["Sample"].unique()
all_plasmids = summary_df["Plasmid"].unique()

# Create full cross-product (all possible sample-plasmid pairs)
all_pairs = pd.MultiIndex.from_product(
    [all_samples, all_plasmids], names=["Sample", "Plasmid"]
)

# Left join to identify missing verdicts
verdict_pairs = verdict_df.set_index(["Sample", "Plasmid"])["VerdictShort"]
full_data = pd.DataFrame(index=all_pairs).join(verdict_pairs, how="left")
full_data["VerdictShort"] = full_data["VerdictShort"].fillna("No Result")

# Aggregate by plasmid
by_plasmid = (
    full_data.reset_index()
    .groupby(["Plasmid", "VerdictShort"])
    .size()
    .unstack(fill_value=0)
)

# Aggregate by sample
by_sample = (
    full_data.reset_index()
    .groupby(["Sample", "VerdictShort"])
    .size()
    .unstack(fill_value=0)
)
```

### Building Hover Data with Name Lists
```python
# Source: Plotly customdata pattern
import pandas as pd

def build_hover_data(verdict_df, group_by_col, list_col):
    """Build dict of verdict -> group -> list of names for hover.

    Args:
        verdict_df: DataFrame with columns [Sample, Plasmid, VerdictShort]
        group_by_col: "Plasmid" or "Sample" (the bar x-axis)
        list_col: "Sample" or "Plasmid" (what to list in hover)

    Returns:
        dict: {verdict: {group_name: [list_of_names]}}
    """
    hover_data = {}
    for verdict in verdict_df["VerdictShort"].unique():
        hover_data[verdict] = {}
        subset = verdict_df[verdict_df["VerdictShort"] == verdict]
        for group_name in subset[group_by_col].unique():
            names = subset[subset[group_by_col] == group_name][list_col].tolist()
            hover_data[verdict][group_name] = names
    return hover_data

# Example usage: grouped by plasmid, list samples
hover_data = build_hover_data(verdict_df, group_by_col="Plasmid", list_col="Sample")
# Result: {"Contaminated": {"pUC19": ["sample1", "sample2"], ...}, ...}
```

### Formatting Hover Text with Truncation
```python
# Source: Plotly hovertemplate + custom logic
def format_hover_text(verdict, count, name_list, max_display=10):
    """Format hover text with verdict, count, and capped name list."""
    if len(name_list) <= max_display:
        names_str = "<br>".join(name_list)
    else:
        names_str = "<br>".join(name_list[:max_display])
        names_str += f"<br>...and {len(name_list) - max_display} more"

    return f"<b>{verdict}</b><br>Count: {count}<br>{names_str}"

# Use in Plotly trace
hover_texts = [
    format_hover_text(verdict, counts[bar], hover_data[verdict].get(bar, []))
    for bar in aggregated_df.index
]

fig.add_trace(go.Bar(
    x=aggregated_df.index,
    y=counts,
    hovertemplate="%{customdata}<extra></extra>",
    customdata=hover_texts,
))
```

### Complete Matplotlib Fallback with Colors
```python
# Source: matplotlib_backend.py pattern + matplotlib stacked bar docs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from plasmicheck.scripts.plotting.colors import HEATMAP_COLORS

def generate_verdict_bar_matplotlib(
    aggregated_df: pd.DataFrame,
    output_path: str,
    title: str,
    xlabel: str,
) -> None:
    """Generate static PNG of verdict stacked bar chart.

    Args:
        aggregated_df: DataFrame with verdict categories as columns, bars as index
        output_path: Output PNG file path
        title: Chart title
        xlabel: X-axis label ("Plasmid" or "Sample")
    """
    if aggregated_df.empty:
        # Handle empty case (same pattern as existing matplotlib functions)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    # Define verdict order and colors (matching Plotly)
    verdict_categories = ["Clean", "Unclear", "Contaminated", "No Result"]
    verdict_colors = {
        "Clean": HEATMAP_COLORS["not_contaminated"],      # lightblue
        "Unclear": HEATMAP_COLORS["unclear"],              # orange
        "Contaminated": HEATMAP_COLORS["contaminated"],    # red
        "No Result": "#D3D3D3",                             # gray
    }

    # Filter to only categories present in data
    present_categories = [c for c in verdict_categories if c in aggregated_df.columns]

    fig, ax = plt.subplots(figsize=(10, 6))

    # Initialize bottom at zero
    bottom = np.zeros(len(aggregated_df))

    # Stack each verdict category
    for verdict in present_categories:
        counts = aggregated_df[verdict].values
        ax.bar(
            aggregated_df.index,
            counts,
            bottom=bottom,
            label=verdict,
            color=verdict_colors[verdict],
            width=0.6,
        )
        bottom += counts

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.legend()

    # Rotate x-axis labels if many items
    if len(aggregated_df) > 5:
        plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Single verdict table | Verdict table + bar charts | Phase 16 (this phase) | Visual overview before diving into tables |
| No missing data handling | Explicit "No Result" category | Phase 16 | Users can distinguish missing from clean |
| plotly.express only | plotly.graph_objects for complex hover | Plotly 5.x best practice | Fine-grained control over customdata |
| Fixed-width charts | Conditional width/rotation for many bars | Current best practice | Readable labels for all dataset sizes |

**Deprecated/outdated:**
None — this is a new feature, no legacy code to replace.

## Open Questions

Things that couldn't be fully resolved:

1. **Scrollbar for 30+ bars**
   - What we know: CSS overflow is standard approach, but keeps axes visible issue
   - What's unclear: Whether to use fixed height + overflow:auto, or just increase chart width
   - Recommendation: Start with increasing chart width proportionally (no scrollbar complexity), revisit if users complain about page width

2. **Optimal "No Result" color shade**
   - What we know: Should be gray to distinguish from verdict colors
   - What's unclear: Exact shade (lightgray #D3D3D3 vs. darkgray #A9A9A9)
   - Recommendation: Use #D3D3D3 (lightgray) matching Backbone_Only color in colors.py

3. **Side-by-side layout on mobile**
   - What we know: Bootstrap col-md-6 used for desktop side-by-side
   - What's unclear: Whether to keep side-by-side on mobile or stack vertically
   - Recommendation: Use Bootstrap responsive grid (col-md-6) — auto-stacks on small screens

## Sources

### Primary (HIGH confidence)
- [Plotly Bar Charts Documentation](https://plotly.com/python/bar-charts/) - Stacked bars with barmode parameter
- [Plotly Hover Text and Formatting](https://plotly.com/python/hover-text-and-formatting/) - customdata and hovertemplate
- [Matplotlib Stacked Bar Chart](https://matplotlib.org/stable/gallery/lines_bars_and_markers/bar_stacked.html) - bottom parameter pattern
- PlasmiCheck codebase (generate_summary_reports.py, matplotlib_backend.py, colors.py) - Existing patterns

### Secondary (MEDIUM confidence)
- [Pandas Stacked Bar Chart Tutorial](https://www.statology.org/pandas-stacked-bar-chart/) - groupby + unstack pattern
- [Python Charts - Stacked Bar Charts](https://www.pythoncharts.com/python/stacked-bar-charts/) - Data preparation
- [Plotly Community - Hovering text in Stacked BarChart](https://community.plotly.com/t/hovering-text-in-stacked-barchart/34813) - Community patterns
- [Plotly Community - Scrollable bar chart](https://community.plotly.com/t/scrollable-bar-chart-or-workarounds/47204) - Handling many bars

### Tertiary (LOW confidence)
- Various Medium/GeeksforGeeks tutorials — provide examples but not authoritative for API details

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — All libraries already in use, verified versions
- Architecture: HIGH — Exact patterns from Phase 13 (metric heatmaps) apply directly
- Pitfalls: MEDIUM-HIGH — Some from official docs (hover, bottom), some inferred (No Result category)

**Research date:** 2026-02-16
**Valid until:** ~90 days (stable libraries, matplotlib/plotly change slowly)
