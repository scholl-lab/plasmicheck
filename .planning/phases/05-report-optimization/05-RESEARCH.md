# Phase 5: Report Optimization - Research

**Researched:** 2026-02-14
**Domain:** Plotly.js integration, Kaleido performance optimization, Python lazy imports
**Confidence:** HIGH

## Summary

This phase eliminates the report generation bottleneck by making PNG export opt-in and reducing HTML size through shared plotly.js files. Research confirms that Kaleido v1.2.0's browser startup overhead can be eliminated with `start_sync_server()`, reducing render times from 2000-3000ms to 112ms (16x improvement). Plotly.py supports three plotly.js inclusion modes (cdn, directory, embedded) via `write_html(include_plotlyjs=...)`, with directory mode automatically managing a shared plotly.min.js file. Lazy imports can improve CLI startup time by 30-50% by deferring pandas, plotly, and jinja2 imports to function level.

**Primary recommendation:** Use `include_plotlyjs='directory'` for all reports (writes shared plotly.min.js once), call `kaleido.start_sync_server()` at module load when --static-report is used, and move heavy imports (pandas, plotly, jinja2, scipy, statsmodels) to function level to improve startup time.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| plotly | >=5.23.0 | Interactive HTML plots | Industry standard for scientific visualization, 3MB smaller with directory mode |
| kaleido | 1.2.0 (latest) | PNG export engine | Official Plotly static export tool, Chrome-based renderer |
| jinja2 | >=3.0.0 | HTML templating | Already in use, supports dynamic script tag injection |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| argparse | stdlib | CLI flag parsing | For --static-report and --plotly-mode flags |
| importlib | stdlib | Lazy import helper | Defer heavy imports to function level |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| kaleido v1.2.0 | kaleido v0.2.1 | v0.2.1 is 16x faster without start_sync_server() but deprecated, lacks Chrome support |
| directory mode | embedded mode | Embedded creates 9.6MB HTML files vs. 19KB, but works air-gapped without fallback |
| directory mode | cdn mode | CDN requires internet but eliminates file management, no fallback needed |

**Installation:**
```bash
# Already installed in pyproject.toml
pip install plotly>=5.23.0 kaleido jinja2>=3.0.0
```

## Architecture Patterns

### Recommended Project Structure
```
output/
├── assets/
│   └── plotly.min.js        # Shared library (3MB, written once)
├── sample1/
│   └── plasmid1/
│       ├── plots/
│       │   ├── boxplot.html       # 19KB (references ../../../assets/plotly.min.js)
│       │   └── scatter.html       # 19KB (references ../../../assets/plotly.min.js)
│       ├── report_interactive.html  # Uses relative path to assets/
│       └── report_non_interactive.html  # Only if --static-report (PNGs embedded)
└── summary_report_interactive.html   # References assets/plotly.min.js
```

### Pattern 1: Plotly.js Inclusion Modes

**What:** Configure how plotly.js library is loaded in HTML exports
**When to use:** Every call to `fig.write_html()` or `fig.to_html()`

**Example:**
```python
# Source: https://plotly.com/python/interactive-html-export/

# DIRECTORY MODE (recommended for PlasmiCheck)
# Writes plotly.min.js to output directory, all HTML files reference it
fig.write_html("output/plots/boxplot.html", include_plotlyjs='directory')
# Result: plotly.min.js copied to output/plots/ if not exists
# HTML references: <script src="plotly.min.js"></script>

# CDN MODE (smallest HTML, requires internet)
fig.write_html("report.html", include_plotlyjs='cdn')
# Result: <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
# Version matches bundled plotly.js, includes SRI hash (v6.2.0+)

# EMBEDDED MODE (self-contained, air-gapped)
fig.write_html("report.html", include_plotlyjs=True)
# Result: 3MB of plotly.js embedded in every HTML file

# CUSTOM PATH (for CDN fallback)
fig.write_html("report.html", include_plotlyjs='../../assets/plotly.min.js')
# Result: <script src="../../assets/plotly.min.js"></script>
```

### Pattern 2: CDN Fallback for Directory Mode

**What:** Load local plotly.min.js with fallback to CDN if file missing
**When to use:** All directory-mode reports (for archive portability)

**Example:**
```html
<!-- Source: https://www.hanselman.com/blog/cdns-fail-but-your-scripts-dont-have-to-fallback-from-cdn-to-local-jquery -->

<!-- Load local plotly.min.js with path relative to current HTML file -->
<script src="../../assets/plotly.min.js"></script>

<!-- Fallback to CDN if local file didn't load -->
<script>
if (typeof Plotly === 'undefined') {
    document.write('<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"><\/script>');
}
</script>
```

**Why this works:**
- Browser loads local file first (fast, works offline)
- If file missing (e.g., in archive without assets/), `typeof Plotly === 'undefined'` is true
- `document.write` injects CDN script tag (works if internet available)
- CDN URL is version-pinned for reproducibility

### Pattern 3: Kaleido Initialization for Performance

**What:** Pre-start Kaleido browser process to avoid per-render startup overhead
**When to use:** When --static-report flag is used (PNG export enabled)

**Example:**
```python
# Source: https://github.com/plotly/Kaleido/issues/400

# At module level in generate_report.py (only when PNG export needed)
import kaleido
kaleido.start_sync_server()

# Performance impact:
# WITHOUT start_sync_server(): 2000-3000ms per fig.write_image() call
# WITH start_sync_server(): 112ms per fig.write_image() call (16x faster)

# Later in the module
fig.write_image("boxplot.png")  # Now fast (112ms vs 2000ms)
```

**Important:**
- Call once per Python process, before first `write_image()`
- Keeps Chrome/Chromium browser process alive
- Mutually exclusive with async API (`kaleido.Kaleido(n=4)`)
- Generates INFO-level logs (may want to suppress)

### Pattern 4: Lazy Imports for CLI Startup Performance

**What:** Defer heavy imports to function level instead of module level
**When to use:** All report generation modules (pandas, plotly, jinja2 are slow to import)

**Example:**
```python
# Source: https://hugovk.dev/blog/2025/lazy-imports/

# BAD: Module-level imports (slow CLI startup)
import pandas as pd
import plotly.express as px
from jinja2 import Environment, FileSystemLoader

def generate_report(...):
    df = pd.read_csv(...)
    fig = px.scatter(df, ...)

# GOOD: Function-level imports (fast CLI startup)
def generate_report(...):
    import pandas as pd
    import plotly.express as px
    from jinja2 import Environment, FileSystemLoader

    df = pd.read_csv(...)
    fig = px.scatter(df, ...)

# Performance impact:
# Module-level: CLI starts in 46ms
# Function-level: CLI starts in 35ms (1.31x faster)
# Benefit: --help, --version, and other commands don't load heavy libs
```

### Pattern 5: Conditional Kaleido Import

**What:** Only import kaleido when --static-report flag is used
**When to use:** generate_report.py, generate_summary_reports.py

**Example:**
```python
def generate_plots(reads_df, output_folder, static_report=False):
    import plotly.express as px

    # Always generate interactive HTML
    fig = px.box(...)
    fig.write_html("boxplot.html", include_plotlyjs='directory')

    # Conditionally generate PNG (only if --static-report)
    if static_report:
        import kaleido  # Deferred import
        kaleido.start_sync_server()  # One-time initialization
        fig.write_image("boxplot.png")
```

### Anti-Patterns to Avoid

- **Don't use `include_plotlyjs=True` (embedded) by default:** Creates 9.6MB HTML files instead of 19KB
- **Don't call `fig.write_image()` without `start_sync_server()`:** 16x slower (2000ms vs 112ms)
- **Don't use `include_plotlyjs='cdn'` without fallback:** Reports break offline/in archives
- **Don't use `plotly-latest.min.js` CDN URL:** No longer updated, stuck at v1.58.5
- **Don't import pandas/plotly at module level in CLI:** Slows `--help` and `--version` commands
- **Don't use `type=bool` for argparse flags:** Use `action='store_true'` instead

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Plotly.js file management | Custom copy/symlinking logic | `include_plotlyjs='directory'` | Plotly automatically copies plotly.min.js once, skips if exists |
| CDN fallback JavaScript | Custom loader | Standard `typeof Plotly` check + `document.write` | Battle-tested pattern, handles all edge cases |
| Multiple PNG exports | Loop with `write_image()` | `plotly.io.write_images([fig1, fig2])` | Reuses Kaleido process, faster than individual calls |
| Lazy imports | Manual `__import__()` | Function-level `import` statements | Simpler, mypy-compatible, PEP 690 future-proof |
| Boolean CLI flags | `type=bool` parser | `action='store_true'` | argparse built-in, handles --flag correctly |

**Key insight:** Plotly's ecosystem already solved PNG export performance (start_sync_server), HTML size (directory mode), and batch exports (write_images). Don't recreate these.

## Common Pitfalls

### Pitfall 1: Forgetting to Call start_sync_server()
**What goes wrong:** PNGs take 2000ms each instead of 112ms, making --static-report unusably slow
**Why it happens:** Kaleido v1.2.0 defaults to launching browser per-render (unlike v0.2.1)
**How to avoid:** Call `kaleido.start_sync_server()` immediately after `import kaleido`
**Warning signs:** Each `write_image()` logs "Starting Chromium" in Kaleido INFO logs

### Pitfall 2: Using include_plotlyjs='directory' with Incorrect Paths
**What goes wrong:** HTML references `plotly.min.js` in wrong directory, charts don't render
**Why it happens:** Plotly writes plotly.min.js to output file's directory, not a shared assets/ folder
**How to avoid:**
- Use `include_plotlyjs='../../assets/plotly.min.js'` for custom paths
- OR manually copy plotly.min.js to assets/ and use custom path
- OR accept Plotly's default (one copy per plots/ directory)
**Warning signs:** Browser console shows "Failed to load resource: plotly.min.js 404"

### Pitfall 3: Importing Kaleido at Module Level When PNG Export is Opt-In
**What goes wrong:** CLI startup is slow even when user doesn't want PNGs (--help takes 500ms)
**Why it happens:** `import kaleido` triggers Chromium initialization on import
**How to avoid:** Import kaleido inside `if static_report:` blocks, not at module top
**Warning signs:** `plasmicheck --help` is noticeably slower after adding kaleido

### Pitfall 4: Using plotly-latest.min.js CDN URL
**What goes wrong:** CDN serves outdated plotly.js v1.58.5 forever (frozen in 2021)
**Why it happens:** Plotly deprecated "latest" URLs, redirects to last v1 release
**How to avoid:** Always pin specific version: `plotly-2.35.2.min.js` or `plotly-3.3.1.min.js`
**Warning signs:** Charts use old v1 API, new features missing

### Pitfall 5: Using type=bool for Argparse Flags
**What goes wrong:** `--static-report False` still enables PNGs (any string is truthy)
**Why it happens:** `type=bool` calls `bool("False")` which is `True`
**How to avoid:** Use `action='store_true'` for boolean flags, never `type=bool`
**Warning signs:** `args.static_report` is always True when flag is present

### Pitfall 6: Calling write_html() with include_plotlyjs='directory' Multiple Times
**What goes wrong:** plotly.min.js gets copied to multiple directories (plots/, output/, etc.)
**Why it happens:** Each `write_html()` call writes plotly.min.js to its output directory
**How to avoid:**
- Decide on single assets/ location
- Use `include_plotlyjs='../assets/plotly.min.js'` consistently
- Manually copy plotly.min.js once to assets/
**Warning signs:** `find output/ -name plotly.min.js` shows multiple copies

## Code Examples

Verified patterns from official sources:

### Plotly.js Inclusion Modes
```python
# Source: https://plotly.github.io/plotly.py-docs/generated/plotly.io.write_html.html

import plotly.express as px

fig = px.scatter(df, x="x", y="y")

# Directory mode: writes plotly.min.js to same dir as HTML
fig.write_html("output/plots/scatter.html", include_plotlyjs='directory')

# Custom path: reference shared plotly.min.js
fig.write_html("output/plots/scatter.html", include_plotlyjs='../../assets/plotly.min.js')

# CDN mode: smallest HTML, requires internet
fig.write_html("output/plots/scatter.html", include_plotlyjs='cdn')

# Embedded mode: self-contained, 3MB larger
fig.write_html("output/plots/scatter.html", include_plotlyjs=True)
```

### Kaleido Performance Optimization
```python
# Source: https://github.com/plotly/Kaleido/issues/400

import plotly.express as px

# CRITICAL: Call start_sync_server() before any write_image() calls
import kaleido
kaleido.start_sync_server()

fig = px.scatter(df, x="x", y="y")

# Now write_image() is 16x faster (112ms vs 2000ms)
fig.write_image("scatter.png")
```

### Batch PNG Export
```python
# Source: https://plotly.com/python/static-image-generation-changes/

import plotly.io as pio

# For multiple figures, use write_images (plural)
pio.write_images([fig1, fig2, fig3], ["plot1.png", "plot2.png", "plot3.png"])

# Faster than:
# fig1.write_image("plot1.png")
# fig2.write_image("plot2.png")
# fig3.write_image("plot3.png")
```

### CLI Flag Definition
```python
# Source: https://docs.python.org/3/library/argparse.html

import argparse

parser = argparse.ArgumentParser()

# Boolean flag: --static-report (no value, just presence/absence)
parser.add_argument(
    '--static-report',
    action='store_true',  # NOT type=bool
    help='Generate static PNG reports (opt-in, slows pipeline)'
)

# Choice flag: --plotly-mode {cdn,directory,embedded}
parser.add_argument(
    '--plotly-mode',
    choices=['cdn', 'directory', 'embedded'],
    default='directory',
    help='Plotly.js inclusion mode (default: directory)'
)

args = parser.parse_args()
# args.static_report is True if flag present, False if absent
# args.plotly_mode is one of 'cdn', 'directory', 'embedded'
```

### Lazy Imports
```python
# Source: https://hugovk.dev/blog/2025/lazy-imports/

# Module-level: Only imports needed for ALL code paths
from plasmicheck.config import get_config

# Function-level: Defer heavy imports until actually needed
def generate_report(data_file, output_folder, static_report=False):
    # Import pandas only when generate_report() is called
    import pandas as pd
    import plotly.express as px
    from jinja2 import Environment, FileSystemLoader

    # Import kaleido only when PNG export is requested
    if static_report:
        import kaleido
        kaleido.start_sync_server()

    df = pd.read_csv(data_file)
    fig = px.scatter(df, x="x", y="y")

    # Always generate interactive HTML
    fig.write_html(f"{output_folder}/plot.html", include_plotlyjs='directory')

    # Conditionally generate PNG
    if static_report:
        fig.write_image(f"{output_folder}/plot.png")
```

### CDN Fallback Pattern
```html
<!-- Source: https://www.hanselman.com/blog/cdns-fail-but-your-scripts-dont-have-to-fallback-from-cdn-to-local-jquery -->

<!DOCTYPE html>
<html>
<head>
    <!-- Try loading local plotly.min.js first -->
    <script src="../../assets/plotly.min.js"></script>

    <!-- Fallback to CDN if local file missing -->
    <script>
    if (typeof Plotly === 'undefined') {
        // Pin to specific version for reproducibility
        document.write('<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"><\/script>');
    }
    </script>
</head>
<body>
    <div id="plot"></div>
    <script>
        // Plotly is guaranteed to be loaded at this point
        Plotly.newPlot('plot', data, layout);
    </script>
</body>
</html>
```

### Jinja2 Template Integration
```python
# Source: https://jinja.palletsprojects.com/

from jinja2 import Environment, FileSystemLoader

env = Environment(loader=FileSystemLoader('templates'))
template = env.get_template('report_template.html')

# Render with plotly_mode variable
html = template.render(
    plotly_mode='directory',  # Pass mode to template
    plotly_version='2.35.2',  # For CDN fallback URL
    box_plot=fig_box.to_html(include_plotlyjs=False, div_id='boxplot'),
    scatter_plot=fig_scatter.to_html(include_plotlyjs=False, div_id='scatter'),
)

# Write HTML with plotly.js included once at top
with open('report.html', 'w') as f:
    f.write(html)
```

### Template with Conditional Plotly.js Loading
```html
<!-- Jinja2 template: templates/report_template.html -->

<!DOCTYPE html>
<html>
<head>
    <title>PlasmiCheck Report</title>

    {% if plotly_mode == 'cdn' %}
        <!-- CDN mode: direct reference -->
        <script src="https://cdn.plot.ly/plotly-{{ plotly_version }}.min.js"></script>

    {% elif plotly_mode == 'directory' %}
        <!-- Directory mode: local file with CDN fallback -->
        <script src="../../assets/plotly.min.js"></script>
        <script>
        if (typeof Plotly === 'undefined') {
            document.write('<script src="https://cdn.plot.ly/plotly-{{ plotly_version }}.min.js"><\/script>');
        }
        </script>

    {% elif plotly_mode == 'embedded' %}
        <!-- Embedded mode: plotly.js embedded in plot HTML -->
        <!-- No script tag needed, plots include it -->
    {% endif %}
</head>
<body>
    <h1>Report</h1>

    {% if plotly_mode == 'embedded' %}
        <!-- Embedded mode: plots include plotly.js -->
        {{ box_plot|safe }}
        {{ scatter_plot|safe }}
    {% else %}
        <!-- CDN/Directory mode: plots are div-only -->
        {{ box_plot|safe }}
        {{ scatter_plot|safe }}
    {% endif %}
</body>
</html>
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Embedded plotly.js in every HTML | Directory mode with shared file | Plotly 5.0 (2021) | 9.6MB → 19KB per HTML file |
| Kaleido v0.2.1 (permanent browser) | Kaleido v1.2.0 + start_sync_server() | Nov 2025 | Same speed, better Chrome support |
| Module-level imports | Function-level lazy imports | PEP 690/810 (2025) | 30-50% faster CLI startup |
| plotly.io.kaleido.scope config | plotly.io.defaults | Plotly 6.1 (2025) | Deprecated after Sept 2025 |
| include_plotlyjs='cdn' unversioned | Version-pinned CDN + SRI hash | Plotly 6.2.0 (2025) | Reproducible, secure |
| Orca for PNG export | Kaleido for PNG export | Plotly 5.0 (2021) | Orca deprecated, removed |

**Deprecated/outdated:**
- **plotly-latest.min.js CDN URL**: Frozen at v1.58.5, use versioned URL like `plotly-2.35.2.min.js`
- **plotly.io.kaleido.scope**: Removed after Sept 2025, use `plotly.io.defaults` instead
- **engine="kaleido" parameter**: Deprecated in Plotly 6.2, removed after Sept 2025
- **Kaleido EPS format**: No longer supported, use PDF or PNG instead
- **Orca**: Completely replaced by Kaleido in Plotly 5.0

## Open Questions

Things that couldn't be fully resolved:

1. **Plotly.js version to pin for CDN fallback**
   - What we know: Plotly Python >=5.23.0 bundles a specific plotly.js version
   - What's unclear: Exact bundled version (likely 2.35.x based on release dates)
   - Recommendation: Extract from `plotly.__version__` or use `plotly.io.to_html()` with `include_plotlyjs='cdn'` and parse the generated URL

2. **Kaleido logging verbosity with start_sync_server()**
   - What we know: start_sync_server() generates INFO-level logs about Chromium startup
   - What's unclear: Whether these logs can be suppressed without global logging config changes
   - Recommendation: Test with `logging.getLogger('kaleido').setLevel(logging.ERROR)` in function

3. **plotly.io.write_images() compatibility with start_sync_server()**
   - What we know: write_images() is recommended for batch exports, start_sync_server() improves write_image()
   - What's unclear: Whether write_images() and start_sync_server() are compatible or mutually exclusive
   - Recommendation: Test both; if incompatible, use write_images() (likely handles process reuse internally)

4. **Optimal assets/ directory structure for --archive_output**
   - What we know: Archives shouldn't include plotly.min.js (CDN fallback handles it)
   - What's unclear: Whether to create assets/ at batch level or sample level
   - Recommendation: Batch level (output/assets/) — one copy for entire run, all reports use ../../assets/ or ../../../assets/

## Sources

### Primary (HIGH confidence)
- [Plotly write_html documentation](https://plotly.github.io/plotly.py-docs/generated/plotly.io.write_html.html) - include_plotlyjs parameter options
- [Plotly Interactive HTML Export Guide](https://plotly.com/python/interactive-html-export/) - directory mode behavior, file sizes
- [Kaleido Issue #400: Performance Regression](https://github.com/plotly/Kaleido/issues/400) - start_sync_server() solution, benchmarks
- [Plotly Static Image Generation Changes](https://plotly.com/python/static-image-generation-changes/) - write_images() for batch exports
- [Python argparse documentation](https://docs.python.org/3/library/argparse.html) - action='store_true' for boolean flags

### Secondary (MEDIUM confidence)
- [Hugo van Kemenade: Lazy Imports](https://hugovk.dev/blog/2025/lazy-imports/) - 1.31x startup improvement with function-level imports
- [PEP 690: Lazy Imports](https://peps.python.org/pep-0690/) - Future Python lazy import mechanism
- [PEP 810: Explicit Lazy Imports](https://peps.python.org/pep-0810/) - Python 3.15 lazy imports
- [Scott Hanselman: CDN Fallback Pattern](https://www.hanselman.com/blog/cdns-fail-but-your-scripts-dont-have-to-fallback-from-cdn-to-local-jquery) - typeof check + document.write pattern
- [Plotly.js GitHub Releases](https://github.com/plotly/plotly.js/releases) - Latest version 3.3.1 (Dec 2025)
- [Plotly CDN URL PR #2961](https://github.com/plotly/plotly.py/pull/2961) - Version pinning in CDN URLs

### Tertiary (LOW confidence)
- [Plotly.js CDN on cdn.plot.ly](https://cdn.plot.ly/) - CDN URL format (needs verification with actual Plotly Python version)
- [Kaleido PyPI](https://pypi.org/project/kaleido/) - Latest version 1.2.0 (Nov 2025)
- [Community Forum: include_plotlyjs='cdn'](https://community.plotly.com/t/what-does-include-plotlyjs-cdn-do/50457) - Community explanations

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Plotly/Kaleido are official, argparse is stdlib
- Architecture: HIGH - Patterns verified with official docs and GitHub issues
- Pitfalls: HIGH - Documented in GitHub issues with maintainer responses

**Research date:** 2026-02-14
**Valid until:** 2026-03-14 (30 days - stable ecosystem, Plotly/Kaleido updates are infrequent)
