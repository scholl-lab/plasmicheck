# PlasmiCheck Report UI/UX Audit

**Auditor:** Senior UI/UX Designer Assessment
**Date:** 2026-02-14
**Scope:** Single-sample HTML report (`report_interactive.html` and `report_template.html`)
**Tested viewports:** 1920x1080 (desktop), 1366x768 (laptop), 768x1024 (tablet)
**Tool:** Playwright headless Chromium screenshots

---

## Overall Ratings

| Aspect | New Report (v0.32.0) | Old Report (v0.25.2) |
|--------|---------------------|---------------------|
| **Visual Design** | 3/10 | 3/10 |
| **Information Architecture** | 5/10 | 4/10 |
| **Responsive Design** | 3/10 | 2/10 |
| **Data Presentation** | 5/10 | 3/10 |
| **Accessibility** | 3/10 | 2/10 |
| **Professional Polish** | 3/10 | 3/10 |
| **Overall** | **3.5/10** | **2.8/10** |

Both reports are "developer-built" — functional but with minimal design thinking. The new report improves data presentation (5-category table, percentages, badges) but introduces layout problems and doesn't address the fundamental design weaknesses.

---

## Critical Issues (Must Fix)

### 1. Raw Python Dict Rendered in HTML
**Severity:** Critical
The `MismatchesNearINSERT` value renders as:
```
{'with_mismatches_or_clipping': 9307, 'without_mismatches_or_clipping': 77}
```
This is a raw Python dictionary dumped into the template. Completely unacceptable for a user-facing report. Should be formatted as two separate labeled values.

### 2. Fixed-Width Plots Don't Scale
**Severity:** Critical
Plots are hard-coded to `width: 1000px` with `max-width: 100%`. On a 1920px display, they sit in the left half with massive empty space. On tablet (768px), they're squeezed. Plotly charts should use `responsive: true` config and be set to `width: 100%` of their container.

### 3. No Visual Hierarchy for Verdict
**Severity:** High
The verdict (the single most important piece of information) is just centered colored text. It should be a prominent banner or card with strong visual contrast. A scientist opening this report wants to know the answer in 1 second.

### 4. Massive White Space on Desktop
**Severity:** High
Content is narrowly centered in a column. The table is `width: 60%`, the `.container` has no max-width. On 1920px, there's ~400px of empty space on each side. The layout wastes space instead of using it purposefully.

---

## Major Issues

### 5. Bootstrap 4.5.2 is Outdated
Bootstrap 5 has been stable since 2021. Bootstrap 4.x has end-of-life concerns. Modernizing to Bootstrap 5 brings better utilities, no jQuery dependency, and improved responsive grid.

### 6. Downsample Notice in Red
Red text signals errors or danger. A downsampling notice is informational. Should use a blue/yellow info banner instead.

### 7. "Plots" is Not a Meaningful Section Header
"Plots" tells the user nothing. Better: "Score Distribution" and "Score Comparison" as individual headings for each chart. Each plot should have a brief explanatory subtitle.

### 8. No Section Navigation
The report requires scrolling through the entire page. For longer reports or when revisiting, there's no way to jump to a specific section. A sticky sidebar or top nav would help.

### 9. Empty Rows Waste Space
When Backbone_Only = 0 and Ambiguous = 0 (which happens frequently), those rows are visual noise. Could collapse or dim them more aggressively.

### 10. No Visual Gauge for Contamination Ratio
The ratio (0.048 vs 4.829) is just a number. A visual gauge, meter, or scale showing where the ratio falls relative to the threshold and unclear range would make the verdict immediately intuitive.

---

## Minor Issues

### 11. Logo Resolution
The base64 PNG logo appears low-resolution on high-DPI displays. Should use SVG for crisp rendering at any scale.

### 12. Metadata Section Lacks Structure
Script version, file paths, and dates are just `<p>` tags in a plain div. These should be in a structured, collapsible details panel.

### 13. No Print Stylesheet
Printing the report would produce poorly formatted output. A `@media print` stylesheet would ensure clean printouts.

### 14. Color Palette Not Colorblind-Safe
The plot colors (default Plotly palette or custom ASSIGNMENT_COLORS) haven't been verified for colorblind accessibility. Scientific reports should use perceptually uniform, colorblind-safe palettes.

### 15. Table Styling Inconsistency
The Read Assignment Summary uses Bootstrap `table-bordered` while the Additional Metrics table uses the same class but looks visually disconnected. They should share a consistent card-based layout.

---

## Responsive Behavior Analysis

| Viewport | Issues |
|----------|--------|
| **1920px (desktop)** | Content sits in narrow center column; massive empty margins; plots fixed 1000px wide |
| **1366px (laptop)** | Slightly better proportions but still wastes space; plots don't fill width |
| **768px (tablet)** | Plots overflow or compress; table columns squeeze; metadata text wraps awkwardly |

---

## Comparison with Industry Standards

### MultiQC (Gold Standard for Bioinformatics Reports)
- Uses collapsible sections with smooth animations
- Sticky navigation sidebar
- Responsive HighCharts plots that fill container width
- Card-based layout with visual grouping
- Consistent, professional color scheme
- Toolbox for customization (hide/show columns, export data)

### FastQC
- Clean pass/warn/fail traffic-light system
- Modular sections with clear visual breaks
- Each section is self-contained with its own plot + summary

### What PlasmiCheck Should Learn
1. **Card-based layout** with clear section boundaries
2. **Responsive plots** that fill their container
3. **Visual verdict indicator** (traffic light, gauge, or prominent banner)
4. **Collapsible sections** for metadata and technical details
5. **Self-contained** (embedded JS) as default

---

## Recommendations (Prioritized)

### Phase 1: Quick Wins (Low Effort, High Impact)
1. Fix MismatchesNearINSERT formatting (split into two labeled values)
2. Make Plotly charts responsive (`responsive: true`, container width 100%)
3. Style verdict as a prominent card/banner with background color
4. Change downsample notice from red to info-blue
5. Rename "Plots" to meaningful section titles

### Phase 2: Layout Overhaul (Medium Effort, High Impact)
1. Upgrade to Bootstrap 5
2. Implement max-width container (1200px) with proper margins
3. Card-based section grouping (verdict card, plots card, data card, metadata card)
4. Add a visual contamination gauge/meter
5. Make metadata collapsible (`<details>` element)

### Phase 3: Polish (Medium Effort, Medium Impact)
1. Add section navigation (sticky top or sidebar)
2. Use colorblind-safe palette
3. Add print stylesheet
4. Convert logo to SVG
5. Add export button (download data as CSV/TSV)

---

## Reference Best Practices

| Practice | Source |
|----------|--------|
| Data tables should combine visual cues (icons, color-coding, badges) | [Pencil & Paper - Data Table UX](https://www.pencilandpaper.io/articles/ux-pattern-analysis-enterprise-data-tables) |
| Use responsive containers, not fixed widths | [CSS-Tricks - Responsive Data Tables](https://css-tricks.com/responsive-data-tables/) |
| Keep visualizations simple with clear message | [TimeTable - Visualization Best Practices](https://www.timetackle.com/data-visualization-best-practices/) |
| Enforce minimum perceptual distance between colors for accessibility | [Accessible Color Sequences (arXiv)](https://arxiv.org/pdf/2107.02270) |
| Scientific reports should use clean layouts, clear hierarchies | [Plotivy - Scientific Visualization Guide](https://plotivy.app/blog/complete-guide-scientific-data-visualization) |
| MultiQC uses cards, collapsible sections, responsive charts | [MultiQC GitHub](https://github.com/MultiQC/MultiQC) |
| Avoid 6+ categorical colors; prefer colorblind-safe palettes | [Plotly Community - Color Selection](https://community.plotly.com/t/how-to-choose-the-right-colors-for-your-data-visualization/86197) |
| Stacking technique for mobile-responsive tables | [NinjaTables - Responsive Tables](https://ninjatables.com/responsive-tables/) |
