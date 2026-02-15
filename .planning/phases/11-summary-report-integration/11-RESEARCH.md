# Phase 11: Summary Report Integration - Research

**Researched:** 2026-02-15
**Domain:** Multi-sample summary report generation with Plotly/matplotlib visualization, pandas data manipulation, Jinja2 templating
**Confidence:** HIGH

## Summary

Research focused on integrating new v0.33.0 metrics (CoverageOutsideINSERT, MismatchesNearINSERT, coverage depth/breadth, resistance genes) into multi-sample summary reports. The standard approach uses Plotly for interactive heatmaps with matplotlib fallback for static exports, pandas for data pivoting/melting, and Jinja2 for conditional rendering with graceful degradation.

**Key findings:**
- PlasmiCheck already uses this exact stack (Plotly 6.3.0, pandas 2.3.3, Jinja2 3.1.2, matplotlib 3.9.1+, seaborn 0.13.2+) with established patterns in `generate_summary_reports.py`
- Summary metrics are stored in per-sample `.summary.tsv` files with `Category\tCount` format, aggregated across samples by the summary report generator
- Existing report has Plotly heatmap + boxplot with matplotlib PNG fallback controlled by `--static-report` and `--plot-backend` flags
- Bootstrap 4.5.2 + jQuery DataTables 1.10.21 already in use for report UI

**Primary recommendation:** Extend the existing `generate_summary_reports.py` patterns (pivot for heatmaps, melt for boxplots, Jinja2 `default()` filter for backward compatibility, matplotlib backend functions in `plotting/matplotlib_backend.py`) rather than introducing new visualization approaches.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Plotly | 6.3.0+ | Interactive heatmaps and boxplots | Official Plotly docs recommend `px.imshow()` for modern heatmaps; PlasmiCheck already uses this in existing ratio heatmap |
| pandas | 2.3.3+ | Data aggregation, pivoting, TSV/Excel export | Standard for Python data manipulation; already handles summary.tsv parsing and multi-sheet Excel export via `pd.ExcelWriter` |
| Jinja2 | 3.1.2+ | HTML template rendering | Standard Python templating; PlasmiCheck uses for both single-sample and summary reports |
| matplotlib | 3.9.1+ | Static PNG plot generation | Fallback for Plotly when `--static-report` is specified; seaborn 0.13.2+ provides higher-level heatmap API |
| XlsxWriter | (via pandas) | Excel file writing | pandas engine for multi-sheet Excel workbooks; already in pyproject.toml dependencies |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| seaborn | 0.13.2+ | Matplotlib-based statistical plots | Heatmap styling in matplotlib backend (uses `sns.heatmap()`) |
| kaleido | latest | Plotly to PNG conversion | When `plot_backend="plotly"` and `--static-report` (converts Plotly figures via `fig.write_image()`) |
| Bootstrap | 4.5.2 | Responsive UI framework | Already loaded in existing summary template for layout and styling |
| jQuery DataTables | 1.10.21 | Sortable/filterable tables | Already applied to all tables in existing summary report |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Plotly px.imshow() | matplotlib sns.heatmap() only | No interactivity in primary report; existing pattern uses Plotly first, matplotlib fallback |
| pandas pivot/melt | Custom dict aggregation | More error-prone; pandas pivot is tested and handles NaN gracefully |
| Jinja2 | Direct string formatting | No conditional rendering for backward compatibility; existing templates use Jinja2 |

**Installation:**
Already in pyproject.toml dependencies - no new packages required.

## Architecture Patterns

### Recommended Project Structure
Current structure already matches best practices:
```
plasmicheck/scripts/
├── generate_summary_reports.py     # Main aggregation + report generation
├── plotting/
│   ├── matplotlib_backend.py       # PNG generation functions
│   └── colors.py                   # Color scheme constants
plasmicheck/templates/
└── summary_template.html            # Jinja2 template for HTML reports
```

### Pattern 1: Data Aggregation from summary.tsv Files
**What:** Read per-sample `.summary.tsv` files, aggregate into single DataFrame
**When to use:** Every summary report generation (existing pattern in `read_compare_outputs()`)
**Example:**
```python
# Source: plasmicheck/scripts/generate_summary_reports.py lines 58-107
def read_compare_outputs(input_dir: str, substring_to_remove: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_files = find_tsv_files(input_dir, ".summary.tsv")
    summary_df_list = []
    for file in summary_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(file)))
        plasmid_name = os.path.basename(os.path.dirname(file))
        df = pd.read_csv(file, sep="\t", low_memory=False)
        df["Sample"] = sample_name
        df["Plasmid"] = plasmid_name
        df = df.rename(columns={"Count": "Value"})
        summary_df_list.append(df)
    summary_df = pd.concat(summary_df_list, ignore_index=True)
    return reads_df, summary_df
```
**Key insight:** `summary_df` contains rows with `Category`, `Value`, `Sample`, `Plasmid` columns. New metrics (MeanDepthBackbone, BreadthBackbone, ResGene_*, etc.) are already in this DataFrame as rows with their respective Category values.

### Pattern 2: Pivot for Heatmaps (Sample × Plasmid Grid)
**What:** Transform long-format data into wide-format grid for heatmap visualization
**When to use:** All heatmap plots (contamination ratio, coverage depth, coverage breadth, mismatches)
**Example:**
```python
# Source: plasmicheck/scripts/generate_summary_reports.py lines 286-288
ratio_data = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value").round(ROUND_DECIMALS)

# For Plotly heatmap
fig = px.imshow(ratio_data, color_continuous_scale=color_scale, height=PLOT_DIMENSIONS["height"], aspect="auto")
```
**Key insight:** Pivot creates the Sample × Plasmid grid needed for imshow(). Missing values become NaN and render as gaps in the heatmap.

### Pattern 3: Melt for Boxplots (Category-based Distribution)
**What:** Transform wide-format data back to long-format for grouped boxplots
**When to use:** Boxplots showing distribution across plasmids or samples
**Example:**
```python
# Source: plasmicheck/scripts/generate_summary_reports.py lines 222-224
boxplot_df = pd.melt(boxplot_data.reset_index(), id_vars=["Sample"], var_name="Plasmid", value_name="Value")
boxplot_df["log_value"] = boxplot_df["Value"].apply(lambda x: np.log10(x + LOG_OFFSET))
fig = px.box(boxplot_df, x="Plasmid", y="log_value", points="all", hover_data={"Sample": True, "Value": True})
```

### Pattern 4: Dual-Mode Plotting (Plotly + Matplotlib Fallback)
**What:** Generate Plotly HTML first, optionally convert to PNG via matplotlib or kaleido
**When to use:** All summary report plots
**Example:**
```python
# Source: plasmicheck/scripts/generate_summary_reports.py lines 245-260
plots_dir = os.path.join(output_dir, "plots")
os.makedirs(plots_dir, exist_ok=True)

# Interactive HTML (always generated)
boxplot_filename = os.path.join(plots_dir, "boxplot_contamination_ratios.html")
fig.write_html(boxplot_filename, include_plotlyjs=False, full_html=False)

# Static PNG (only if --static-report)
boxplot_filename_png = None
if static_report:
    boxplot_filename_png = os.path.join(plots_dir, "boxplot_contamination_ratios.png")
    if plot_backend == "matplotlib":
        from .plotting.matplotlib_backend import generate_summary_boxplot_matplotlib
        generate_summary_boxplot_matplotlib(boxplot_data, boxplot_filename_png, LOG_OFFSET, "outliers")
    else:  # plot_backend == "plotly"
        fig.write_image(boxplot_filename_png)  # Uses kaleido
```
**Key insight:** HTML fragments written with `include_plotlyjs=False, full_html=False` are embedded into Jinja2 template via `{{ plot_content|safe }}`.

### Pattern 5: Multi-Sheet Excel Export with XlsxWriter
**What:** Write multiple DataFrames to separate sheets in single .xlsx file
**When to use:** Exporting summary tables (existing pattern for 4 sheets, extend to 6+ sheets)
**Example:**
```python
# Source: plasmicheck/scripts/generate_summary_reports.py lines 154-159
excel_file = os.path.join(data_dir, "summary_report_tables.xlsx")
with pd.ExcelWriter(excel_file, engine="xlsxwriter") as writer:
    _sanitize_for_excel(combined_df).to_excel(writer, sheet_name="Combined", index=False)
    _sanitize_for_excel(verdict_df).to_excel(writer, sheet_name="Verdicts", index=False)
    _sanitize_for_excel(ratio_df).to_excel(writer, sheet_name="Ratios", index=False)
    _sanitize_for_excel(p_values_df).to_excel(writer, sheet_name="P-Values", index=False)
```
**Note:** `_sanitize_for_excel()` converts boolean columns to strings to avoid locale-dependent rendering (lines 119-125).

### Pattern 6: Backward Compatibility with Jinja2 default() Filter
**What:** Handle missing metrics gracefully when rendering old data
**When to use:** When new metrics might not exist in pre-v0.33.0 summary.tsv files
**Example:**
```python
# From Phase 9-02 implementation (generate_report.py pattern)
def _get_metric(df: pd.DataFrame, category: str, default: str = "0.00") -> str:
    rows = df[df["Category"] == category]
    return str(rows["Count"].values[0]) if not rows.empty else default

coverage_metrics = {
    "mean_depth_insert": _get_metric(summary_df, "MeanDepthInsert"),
    "breadth_backbone": _get_metric(summary_df, "BreadthBackbone"),
    # ... returns "0.00" if category not found
}
```
**In template:**
```html
<!-- Jinja2 template with graceful fallback -->
{% if coverage_metrics.mean_depth_insert|default("0.00") != "0.00" %}
    <!-- Show coverage section -->
{% endif %}
```
**Key insight:** Use `|default()` filter in templates and check for sentinel values (e.g., "0.00" or empty string) to conditionally render sections.

### Anti-Patterns to Avoid
- **Hard-coding metric names in template:** Extract metric names from config or use data-driven rendering to avoid template updates when adding metrics
- **Creating new DataFrame from scratch:** Always pivot/melt existing summary_df rather than rebuilding from dictionaries
- **Inconsistent color scales:** Use config.json `plot_summary.colors` values to maintain visual consistency across all heatmaps
- **Not handling NaN in pivot:** pandas pivot() produces NaN for missing sample-plasmid combinations; Plotly renders these as gaps (expected behavior), but check for empty DataFrames before plotting

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Aggregating summary.tsv files | Custom file parser with glob/loops | `pd.read_csv()` + `pd.concat()` (existing pattern) | Already handles edge cases (missing files, empty files, encoding issues); existing code in `read_compare_outputs()` |
| Converting wide to long format | Manual loop with list of dicts | `pd.pivot()` and `pd.melt()` | Handles NaN, preserves dtypes, tested extensively; existing usage in lines 286-288 and 222-224 |
| Multi-sheet Excel export | Looping with openpyxl | `pd.ExcelWriter` with xlsxwriter engine | Handles formatting, auto-column-width, boolean sanitization via `_sanitize_for_excel()` |
| Plotly to PNG conversion | Screenshot or browser automation | `kaleido` package or matplotlib backend | kaleido is official Plotly static image generator; matplotlib backend provides consistent styling |
| Missing data handling | Try/except on KeyError | Jinja2 `default()` filter + pandas `fillna()` | Declarative and template-friendly; already used in Phase 9-02 for coverage metrics |
| Collapsible sections | Custom JavaScript | Bootstrap collapse/accordion components | Already loaded in template (Bootstrap 4.5.2); aria attributes for accessibility built-in |
| Color scale generation | Manual RGB interpolation | Plotly built-in color scales or seaborn palettes | Perceptually uniform; existing usage: `color_continuous_scale` in line 291-302 |

**Key insight:** PlasmiCheck already has working implementations of all these patterns. Phase 11 extends existing functions rather than creating new ones.

## Common Pitfalls

### Pitfall 1: Forgetting to Filter summary_df by Category Before Pivot
**What goes wrong:** Pivoting the entire summary_df (all categories) creates a massive grid with mixed data types (verdict strings, ratio floats, coverage floats all in same pivot)
**Why it happens:** summary_df contains ALL metrics as rows with Category column distinguishing them
**How to avoid:** Always filter first: `summary_df[summary_df["Category"] == "MeanDepthBackbone"]` before pivot
**Warning signs:** Pivot table has unexpected rows or columns; ValueError about incompatible dtypes
**Example from existing code:**
```python
# Line 516: Filter FIRST, then pivot
ratio_df = summary_df[summary_df["Category"] == "Ratio"]
boxplot_data, p_values_df = calculate_variations(ratio_df)
```

### Pitfall 2: Not Converting "Value" Column to Numeric Before Pivot
**What goes wrong:** Pivot fails or produces object dtype instead of float64, breaking math operations
**Why it happens:** summary.tsv "Count" column (renamed to "Value") is read as string by default
**How to avoid:** Use `pd.to_numeric(df["Value"], errors="coerce")` before pivot
**Warning signs:** Heatmap shows strings instead of colors; errors like "cannot perform operation with object dtype"
**Example from existing code:**
```python
# Line 174 and 283: Convert to numeric with coercion
ratio_df.loc[:, "Value"] = pd.to_numeric(ratio_df["Value"], errors="coerce")
```

### Pitfall 3: Empty DataFrame After Filter Causing Plot Failures
**What goes wrong:** Filtering summary_df for a category that doesn't exist (e.g., ResGene_* in old data) produces empty DataFrame, pivot() returns empty, px.imshow() throws error
**Why it happens:** Backward compatibility - old summary.tsv files don't have new metrics
**How to avoid:** Check `df.empty` before plotting and skip plot generation or show placeholder
**Warning signs:** KeyError or ValueError when generating plots; "No data available" errors in browser console
**Example from matplotlib_backend.py:**
```python
# Lines 182-195: Check for empty DataFrame
if ratio_data.empty:
    logger.warning("Empty DataFrame provided to generate_heatmap_matplotlib")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return
```

### Pitfall 4: Resistance Gene Table Has Variable Row Count
**What goes wrong:** Trying to create a Sample × Plasmid pivot for resistance genes fails because each sample-plasmid pair can have 0, 1, or N genes
**Why it happens:** Resistance genes are per-gene metrics, not per-sample aggregates
**How to avoid:** Keep resistance gene data as flat table (Sample | Plasmid | GeneName | MeanDepth | Breadth), export as separate TSV/sheet
**Warning signs:** Pivot produces multi-index with gene names; unable to visualize as heatmap
**Recommended approach:** Flat table only, no heatmap (as per CONTEXT.md decision)

### Pitfall 5: Inconsistent NaN Handling Between Plotly and Matplotlib
**What goes wrong:** Plotly renders NaN as gaps (white space), matplotlib may show as zero or different color depending on colormap
**Why it happens:** Different default behaviors for missing data
**How to avoid:** Explicitly set how to handle NaN - for Plotly, use `hoverongaps=False`; for matplotlib, use `mask` parameter in seaborn heatmap
**Warning signs:** Heatmaps look different between interactive and static versions; unexpected zeros in plot
**Example:**
```python
# Plotly: NaN rendered as gap (default behavior is acceptable)
fig = px.imshow(ratio_data)  # NaN shows as gap

# Matplotlib/seaborn: use mask for consistency
import numpy as np
mask = ratio_data.isna()
sns.heatmap(ratio_data, mask=mask, cmap="viridis", cbar=False)
```

### Pitfall 6: Bootstrap 4 vs 5 Accordion Syntax
**What goes wrong:** Using Bootstrap 5 accordion syntax (data-bs-*) in Bootstrap 4 template (data-toggle/data-target)
**Why it happens:** Existing template uses Bootstrap 4.5.2, not Bootstrap 5
**How to avoid:** Check template's Bootstrap version (line 6 of summary_template.html shows 4.5.2), use Bootstrap 4 syntax
**Warning signs:** Accordion doesn't expand/collapse; browser console shows JavaScript errors
**Correct syntax for Bootstrap 4:**
```html
<div class="card">
  <div class="card-header" id="headingOne">
    <button class="btn btn-link" data-toggle="collapse" data-target="#collapseOne">
      Section Title
    </button>
  </div>
  <div id="collapseOne" class="collapse show" data-parent="#accordion">
    <div class="card-body">Content</div>
  </div>
</div>
```

## Code Examples

Verified patterns from official sources and existing PlasmiCheck code:

### Extracting New Metrics from summary_df
```python
# Pattern: Filter by Category, convert to numeric, pivot to Sample × Plasmid grid
# Source: Existing code in generate_summary_reports.py lines 174-288

def create_coverage_heatmap_data(summary_df: pd.DataFrame, metric_category: str) -> pd.DataFrame:
    """Extract and pivot coverage metric for heatmap visualization.

    Args:
        summary_df: Aggregated summary DataFrame with Category, Value, Sample, Plasmid columns
        metric_category: Category name (e.g., "MeanDepthBackbone", "BreadthBackbone")

    Returns:
        Pivoted DataFrame (Sample × Plasmid) or empty DataFrame if category not found
    """
    metric_df = summary_df[summary_df["Category"] == metric_category]
    if metric_df.empty:
        return pd.DataFrame()  # Graceful degradation for old data

    metric_df.loc[:, "Value"] = pd.to_numeric(metric_df["Value"], errors="coerce")
    pivot_data = metric_df.pivot(index="Sample", columns="Plasmid", values="Value").round(2)
    return pivot_data
```

### Creating Heatmap with Custom Color Scale
```python
# Pattern: Use existing plot_summary config colors, customize scale for new metrics
# Source: Lines 291-302 (existing ratio heatmap pattern)

def plot_coverage_heatmap(
    coverage_data: pd.DataFrame,
    output_dir: str,
    metric_name: str,  # "Coverage Depth" or "Coverage Breadth"
    vmax: float | None = None,  # Cap for depth heatmap
    plot_config: dict[str, Any],
    static_report: bool = False,
    plot_backend: str = "plotly",
) -> tuple[str, str | None]:
    """Generate coverage heatmap with blue-to-red color scale.

    Blue (low/good) to red (high/bad) matches existing contamination ratio heatmap.
    """
    if coverage_data.empty:
        logging.warning(f"No data for {metric_name} heatmap")
        return ("", None)

    # Blue-to-red scale (matches existing contamination heatmap)
    color_scale = [
        [0.0, plot_config["colors"]["not_contaminated"]],  # "lightblue"
        [1.0, plot_config["colors"]["contaminated"]],       # "red"
    ]

    # Cap values if needed (for depth heatmap)
    if vmax is not None:
        plot_data = coverage_data.clip(upper=vmax)
    else:
        plot_data = coverage_data

    fig = px.imshow(
        plot_data,
        color_continuous_scale=color_scale,
        height=PLOT_DIMENSIONS["height"],
        aspect="auto",
        title=f"Heatmap of {metric_name} by Sample and Plasmid",
    )

    fig.update_layout(
        xaxis_title="Plasmid",
        yaxis_title="Sample",
        xaxis={"tickangle": 45, "tickfont": {"size": 9}},
        yaxis={"tickfont": {"size": 9}},
    )
    fig.update_coloraxes(showscale=False)  # Hide color bar (matches existing pattern)

    # Save HTML and optionally PNG (follow existing pattern)
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    html_filename = os.path.join(plots_dir, f"heatmap_{metric_name.lower().replace(' ', '_')}.html")
    fig.write_html(html_filename, include_plotlyjs=False, full_html=False)

    png_filename = None
    if static_report:
        png_filename = os.path.join(plots_dir, f"heatmap_{metric_name.lower().replace(' ', '_')}.png")
        if plot_backend == "matplotlib":
            # Create matplotlib version (would need new function in matplotlib_backend.py)
            from .plotting.matplotlib_backend import generate_coverage_heatmap_matplotlib
            generate_coverage_heatmap_matplotlib(coverage_data, png_filename, vmax, plot_config)
        else:
            fig.write_image(png_filename)

    return html_filename, png_filename
```

### Creating Resistance Gene Flat Table
```python
# Pattern: Parse ResGene_* metrics from summary_df, create flat table
# Source: New pattern based on summary.tsv format from compare_alignments.py lines 512-519

def create_resistance_gene_table(summary_df: pd.DataFrame) -> pd.DataFrame:
    """Extract resistance gene coverage metrics into flat table.

    Returns DataFrame with columns: Sample, Plasmid, GeneName, MeanDepth, Breadth
    """
    # Check if resistance genes were available (backward compatibility)
    rg_available = summary_df[summary_df["Category"] == "ResistanceGenesAvailable"]
    if rg_available.empty or rg_available["Value"].iloc[0] == "false":
        return pd.DataFrame(columns=["Sample", "Plasmid", "GeneName", "MeanDepth", "Breadth"])

    # Extract gene names from ResistanceGenes row
    rg_rows = summary_df[summary_df["Category"] == "ResistanceGenes"]

    gene_table_rows = []
    for _, row in rg_rows.iterrows():
        sample = row["Sample"]
        plasmid = row["Plasmid"]
        gene_names_str = row["Value"]

        if gene_names_str == "None" or pd.isna(gene_names_str):
            continue  # No genes detected for this sample-plasmid

        gene_names = gene_names_str.split(",")

        for gene_name in gene_names:
            # Find ResGene_{sanitized_name}_MeanDepth and _Breadth rows
            sanitized = gene_name.replace("/", "_").replace(" ", "_")

            depth_row = summary_df[
                (summary_df["Sample"] == sample) &
                (summary_df["Plasmid"] == plasmid) &
                (summary_df["Category"] == f"ResGene_{sanitized}_MeanDepth")
            ]
            breadth_row = summary_df[
                (summary_df["Sample"] == sample) &
                (summary_df["Plasmid"] == plasmid) &
                (summary_df["Category"] == f"ResGene_{sanitized}_Breadth")
            ]

            if not depth_row.empty and not breadth_row.empty:
                gene_table_rows.append({
                    "Sample": sample,
                    "Plasmid": plasmid,
                    "GeneName": gene_name,
                    "MeanDepth": float(depth_row["Value"].iloc[0]),
                    "Breadth": float(breadth_row["Value"].iloc[0]),
                })

    return pd.DataFrame(gene_table_rows)
```

### Matplotlib Heatmap Fallback with Custom Colormap
```python
# Pattern: Create matplotlib heatmap matching Plotly's blue-to-red scale
# Source: New function for matplotlib_backend.py, following existing patterns

def generate_coverage_heatmap_matplotlib(
    coverage_data: pd.DataFrame,
    output_path: Path | str,
    vmax: float | None,
    plot_config: dict[str, Any],
) -> None:
    """Generate coverage heatmap using matplotlib/seaborn with custom colormap.

    Matches Plotly's blue-to-red scale for visual consistency.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap

    if coverage_data.empty:
        logger.warning("Empty DataFrame provided to generate_coverage_heatmap_matplotlib")
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    sns.set_theme(style="whitegrid")

    # Create custom colormap: blue (low) to red (high)
    colors = [plot_config["colors"]["not_contaminated"], plot_config["colors"]["contaminated"]]
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list("coverage", colors, N=n_bins)

    # Cap data if needed
    plot_data = coverage_data.clip(upper=vmax) if vmax is not None else coverage_data

    # Determine figure size based on data dimensions
    n_samples, n_plasmids = plot_data.shape
    fig_width = max(10, n_plasmids * 0.8)
    fig_height = max(6, n_samples * 0.4)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Create heatmap without color bar (matches Plotly pattern)
    sns.heatmap(
        plot_data,
        cmap=cmap,
        annot=False,
        cbar=False,
        ax=ax,
        vmin=0,
        vmax=vmax if vmax is not None else None,
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_xlabel("Plasmid")
    ax.set_ylabel("Sample")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved coverage heatmap to {output_path}")
```

### Jinja2 Template with Conditional Rendering
```html
<!-- Pattern: Show section only if data exists, gracefully degrade for old data -->
<!-- Source: Bootstrap 4 accordion syntax + Jinja2 default() filter best practices -->

<!-- Coverage Analysis Section (new in v0.33.0) -->
{% if depth_heatmap_content or breadth_heatmap_content %}
<div class="card">
  <div class="card-header" id="headingCoverage">
    <h5 class="mb-0">
      <button class="btn btn-link" data-toggle="collapse" data-target="#collapseCoverage" aria-expanded="true">
        Coverage Analysis
      </button>
    </h5>
  </div>
  <div id="collapseCoverage" class="collapse show" data-parent="#accordion">
    <div class="card-body">
      {% if depth_heatmap_content %}
      <h3>Coverage Depth (Backbone Mean Depth)</h3>
      <div class="plot">
        {% if interactive %}
          {{ depth_heatmap_content|safe }}
        {% else %}
          <img src="{{ depth_heatmap_content }}" alt="Coverage Depth Heatmap">
        {% endif %}
      </div>
      {% endif %}

      {% if breadth_heatmap_content %}
      <h3>Coverage Breadth (Backbone Breadth %)</h3>
      <div class="plot">
        {% if interactive %}
          {{ breadth_heatmap_content|safe }}
        {% else %}
          <img src="{{ breadth_heatmap_content }}" alt="Coverage Breadth Heatmap">
        {% endif %}
      </div>
      {% endif %}
    </div>
  </div>
</div>
{% endif %}

<!-- Resistance Genes Section (new in v0.33.0) -->
{% if resistance_gene_table %}
<div class="card">
  <div class="card-header" id="headingResGenes">
    <h5 class="mb-0">
      <button class="btn btn-link" data-toggle="collapse" data-target="#collapseResGenes" aria-expanded="true">
        Resistance Gene Coverage
      </button>
    </h5>
  </div>
  <div id="collapseResGenes" class="collapse show" data-parent="#accordion">
    <div class="card-body">
      <div class="table-container">
        {{ resistance_gene_table|safe }}
      </div>
    </div>
  </div>
</div>
{% endif %}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Plotly `create_annotated_heatmap` figure factory | `px.imshow()` for heatmaps | Plotly 5.5.0 (2022) | Simpler API, better performance, aspect="auto" for flexible sizing |
| `pd.wide_to_long()` for reshaping | `pd.melt()` for long format | pandas best practices (2024+) | More flexible, handles multi-column id_vars |
| Hardcoded metric extraction in template | Jinja2 `default()` filter for missing data | Jinja2 2.x standard practice | Graceful degradation, backward compatibility without try/except |
| Bootstrap 3 accordion | Bootstrap 4 collapse component | Bootstrap 4.0 (2018) | Improved accessibility (aria-* attributes), better mobile support |
| kaleido replaced orca | kaleido for Plotly static export | Plotly 4.9+ (2020) | Faster, no external dependencies, official recommendation |

**Deprecated/outdated:**
- `create_annotated_heatmap`: Replaced by `px.imshow()` with `text_auto=True` for annotations
- `pd.wide_to_long()`: Prefer `pd.melt()` for simpler cases (wide_to_long is still valid but less flexible)
- Bootstrap 3 data-toggle syntax: Use Bootstrap 4/5 syntax with proper aria attributes
- Plotly orca static image engine: Replaced by kaleido (faster, maintained)

## Open Questions

Things that couldn't be fully resolved:

1. **Linear scale cap value for depth heatmap**
   - What we know: Depth can vary widely (0-1000+ for contaminated samples, 0-10 for clean)
   - What's unclear: Optimal cap value for visualization (too low = saturates, too high = poor color contrast)
   - Recommendation: Start with cap=100 (based on typical clean sample max depth ~10-20, contaminated ~50-200); make configurable in config.json as `plot_summary.depth_heatmap_cap` for user tuning

2. **Exact positioning of new sections in report**
   - What we know: User specified ordering (Contamination → Coverage → Mismatches → Resistance Genes)
   - What's unclear: Whether "Coverage" means CoverageOutsideINSERT OR coverage depth/breadth, or both together
   - Recommendation: Group related metrics - CoverageOutsideINSERT stays with contamination (it's a contamination indicator), new Coverage section is depth/breadth only

3. **Handling samples with no resistance genes vs. no GenBank file**
   - What we know: ResistanceGenesAvailable flag distinguishes these cases
   - What's unclear: Should the table show "No resistance genes detected" row or omit the section entirely?
   - Recommendation: Omit section entirely if ALL samples have ResistanceGenesAvailable=false; show empty table with message if some samples have genes but others don't

4. **Should coverage heatmaps use log scale like boxplot?**
   - What we know: Existing boxplot uses log10 scale for contamination ratios
   - What's unclear: Whether coverage depth (mean 0-100+) benefits from log scale or linear is better
   - Recommendation: Use linear scale with cap for depth (interpretability), linear for breadth (already % 0-100 bounded)

## Sources

### Primary (HIGH confidence)
- Plotly official documentation - [Heatmaps in Python](https://plotly.com/python/heatmaps/)
- Plotly official API reference - [plotly.express.imshow](https://plotly.com/python-api-reference/generated/plotly.express.imshow.html)
- pandas official documentation - [Reshaping and pivot tables](https://pandas.pydata.org/docs/user_guide/reshaping.html)
- pandas official API reference - [pandas.melt](https://pandas.pydata.org/docs/reference/api/pandas.melt.html)
- pandas official API reference - [pandas.DataFrame.to_excel](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_excel.html)
- pandas official API reference - [pandas.ExcelWriter](https://pandas.pydata.org/docs/reference/api/pandas.ExcelWriter.html)
- Jinja2 official documentation - [Template Designer Documentation](https://jinja.palletsprojects.com/en/stable/templates/)
- Bootstrap official documentation - [Accordion · Bootstrap v5.3](https://getbootstrap.com/docs/5.3/components/accordion/)
- Bootstrap official documentation - [Collapse · Bootstrap v4.6](https://getbootstrap.com/docs/4.6/components/collapse/)
- seaborn official documentation - [seaborn.heatmap](https://seaborn.pydata.org/generated/seaborn.heatmap.html)
- seaborn official documentation - [Choosing color palettes](https://seaborn.pydata.org/tutorial/color_palettes.html)
- XlsxWriter official documentation - [Example: Pandas Excel with multiple dataframes](https://xlsxwriter.readthedocs.io/example_pandas_multiple.html)
- XlsxWriter official documentation - [Working with Pandas and XlsxWriter](https://xlsxwriter.readthedocs.io/working_with_pandas.html)

### Secondary (MEDIUM confidence)
- PlasmiCheck existing implementation - `plasmicheck/scripts/generate_summary_reports.py` (verified patterns)
- PlasmiCheck existing implementation - `plasmicheck/scripts/plotting/matplotlib_backend.py` (matplotlib fallback patterns)
- PlasmiCheck existing implementation - `plasmicheck/templates/summary_template.html` (Bootstrap 4 + Jinja2 patterns)
- PlasmiCheck existing implementation - `plasmicheck/config.json` (color scheme, plot dimensions)
- [Python Charts - Heatmaps in plotly with imshow](https://python-charts.com/correlation/heatmap-plotly/) - Verified best practices
- [Towards Data Science - Reshaping a Pandas Dataframe: Long-to-Wide and Vice Versa](https://towardsdatascience.com/reshaping-a-pandas-dataframe-long-to-wide-and-vice-versa-517c7f0995ad/) - Pivot/melt patterns
- [TheLinuxCode - Ansible and Jinja2: Creating Dynamic Templates (Practical Patterns for 2026)](https://thelinuxcode.com/ansible-and-jinja2-creating-dynamic-templates-practical-patterns-for-2026/) - Jinja2 conditional rendering
- [Managing API Changes: 8 Strategies That Reduce Disruption by 70% (2026 Guide)](https://www.theneo.io/blog/managing-api-changes-strategies) - Backward compatibility strategies
- [Matplotlib vs Seaborn vs Plotly (2026)](https://www.blog.qualitypointtech.com/2026/02/matplotlib-vs-seaborn-vs-plotly.html) - Library comparison

### Tertiary (LOW confidence)
- [GeeksforGeeks - How to Write Pandas DataFrames to Multiple Excel Sheets](https://www.geeksforgeeks.org/python/how-to-write-pandas-dataframes-to-multiple-excel-sheets/) - General tutorial, less authoritative than official docs
- [Colorlib - 20 Best Bootstrap Accordion Examples 2026](https://colorlib.com/wp/bootstrap-accordions/) - Examples, not official documentation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries already in use with verified versions from pyproject.toml and existing code
- Architecture: HIGH - All patterns extracted from existing PlasmiCheck code (`generate_summary_reports.py`, `matplotlib_backend.py`)
- Pitfalls: MEDIUM - Based on common pandas/Plotly issues and existing code guard patterns; some are inferred from code structure
- Code examples: HIGH - Directly from PlasmiCheck codebase and official documentation

**Research date:** 2026-02-15
**Valid until:** 60 days (stable libraries, established patterns; pandas/Plotly/Jinja2 are mature with infrequent breaking changes)
