# Cohort tables, week labels and timecourse plots -- design

Date: 2026-11-20
Status: approved 2026-11-20

## Intent

The DB-backed analysis docs (first user: `analysis_26-1_db.Rmd`, then the
`analysis-db` template) need three things the package does not yet provide:

1. **A view of the whole mouse cohort.** Today the doc shows `count(group, sex)`
   plus three scattered QC chunks (missing panels, dropouts, exclusions). Janis
   wants (A) a per-group cohort summary that can be quoted in a lab meeting or
   methods section and (B) a per-mouse grid to check that every mouse is
   accounted for across assays.
2. **"Week after BMT" labels** for `collected_at`, so serial samples (titer
   bleeds) and end-points are named and ordered by study week.
3. **Timecourse plots of ELISA results** (anti-MPO titer over time): one with
   the individual mice, one with only the group median and spread.

## Goals

- Four exported functions plus one internal-sharing plot pair:
  `meta_summary_table()`, `meta_mouse_table()`, `meta_label_weeks()`,
  `plot_timecourse()`, `plot_timecourse_summary()`.
- All work on `db_load_experiment()` output (`dat$meta`, `dat$exclusions`, the
  per-unit tibbles) with no DB connection.
- The cohort tables include mice that never reached the DB (see Constraints).
- Wire them into `analysis_26-1_db.Rmd` and the `analysis-db` template.

## Non-goals

- Storing dropouts (mice without `mouse_ID`, early-death reasons) in the DB.
  They stay in the meta workbook, which the doc reads.
- Statistics tables for timecourses (beyond the optional per-week stars).
- Changing `plot_box()` / `plot_population()` / `plot_paired()` behaviour.

## Anti-goals

- A cohort summary that silently reports only the 17 DB mice as "the cohort"
  when 23 were enrolled.
- Week labels that split one end-point harvest across two weeks.
- A summary plot whose spread measure overstates precision (SEM) or draws
  interpolated values between sampling weeks (ribbons).

## Constraints / facts found

- `subjects` (via `db_query_subjects()` -> `dat$meta`) has `experiment_code`,
  `mouse_id`, `group`, `mouse_treatment`, `sex`, `cage`, `mouse_strain`,
  `generation`, `dob`, `start_date`, `bmt_date`, `end_date`. 26-1: all
  `bmt_date` 2026-05-12, end dates 2026-06-10 / 06-15 / 06-16.
- The meta workbook `meta` sheet for 26-1 lists **23** mice; **17** are in the
  DB. The 6 others have no `mouse_ID` (only the `ID` ear tag / MDC id):
  1 died at the blood draw before BMT (`bmt_date` empty, group `irrel`), 2 were
  removed at BMT for a poor titer (`group` empty), 3 died after BMT (05-22,
  05-31, 05-31). One DB mouse (26-1-17, irrel) is an early death at day 29
  (`early = "yes"`, `reason = "Schmerzen"`). Workbook columns: `mouse_ID`,
  `ID`, `cage`, `line`, `generation`, `sex`, `group`, `DOB`, `start_date`,
  `bmt_date`, `death_date`, `early`, `reason`, `comment`.
- `dat$exclusions`: current imports' `import_exclusions` rows (26-1-16
  excluded at export). `db_load_experiment()` drops mice excluded mouse-wide in
  every unit from `dat$meta` unless `keep_excluded = TRUE`.
- Titer data (`anti-mpo-titer_26-1_long`): one replicate, `lloq`/`uloq` `NA`,
  dates 05-05 (d -7), 05-19 (d 7), 06-02 (d 21), end-points 06-10 (d 29, n=1),
  06-15 (d 34), 06-16 (d 35). Restim ELISA: 2 replicates, `OOR<` rows have
  `value = NA` (never substituted -- DB rule).
- Package rules (CLAUDE.md): `meta_` and `plot_` domains exist; `|>`,
  namespaced calls, `.data` only for variable column names, roxygen
  `@param/@returns/@export`, non-ASCII escaped, `devtools::check()` 0/0.
  `gt`, `rstatix`, `ggplot2`, `scales` are already in Imports.

## Approaches considered

Cohort view (agreed with Janis: build both, two functions):

- **A. Per-group summary ("Table 1")** -- characteristics x group. Quotable,
  hides per-mouse detail.
- **B. Per-mouse accounting grid** -- one row per mouse x assay presence and
  status. Replaces three QC chunks; long, internal only.
- One combined function with a `type =` switch was rejected: the inputs and
  outputs differ enough that one signature would be half unused arguments.

Week labels (Janis decided the rule):

- Plain floor/round weeks for every date -- splits the 26-1 end-point
  (d 34 / d 35) only under floor, and any harvest straddling a rounding
  boundary in general.
- "pre" / "end" labels -- rejected by Janis: every measurement gets a week.
- **Chosen:** nearest whole week for every sample; end-point samples from the
  same calendar week are pooled to one week; an end-point in a different week
  (26-1-17, d 29) sits at its own week (week 4).

Timecourse:

- **Chosen:** two functions sharing internals -- `plot_timecourse()` (each
  mouse + group median) and `plot_timecourse_summary()` (median + IQR only).
  A single function with `summary_only = TRUE` was rejected: half its
  arguments would be dead in each mode.
- Spread measure for the summary: **median + IQR** (chosen: skewed titers,
  n ~ 8, matches the rank-based tests). Mean +/- SEM (precision of the mean,
  shrinks with n, mixes frameworks with Wilcoxon stars) and median 95% CI
  (spans almost the whole range at n = 8) rejected.
- Error bars (chosen) vs ribbon: few, unevenly spaced timepoints and two
  overlapping groups -> error bars, groups dodged.

## Design

### 1. `meta_label_weeks(data, meta, by = c("experiment_code", "mouse_id"), ref_col = "bmt_date", date_col = "collected_at")`

File `R/meta_weeks.R`. Adds three columns to `data` and returns it:

- `days_post_bmt` (integer) = `collected_at - bmt_date` of the mouse
  (joined from `meta` by `by`; `meta` = `dat$meta`).
- `week_num` (integer) = `round(days_post_bmt / 7)`, except pooled end-points
  (below). Negative weeks allowed (`-1` = pre-BMT bleed); week 0 possible.
- `week` = ordered factor `"week <n>"`, levels sorted numerically
  (`"week -1" < "week 1" < "week 3" < ...`), only levels that occur.

End-point pooling: a sample is an end-point when its `collected_at` equals the
mouse's `end_date`. Within each `experiment_code`, end-point samples are grouped
by calendar week (ISO week, Mon-Sun, of `collected_at`); each group gets
`week_num = round(median(days_post_bmt) / 7)`. Serial (non-end-point) samples
are never pooled. 26-1 result: d -7 -> -1, d 7 -> 1, d 21 -> 3,
d 29 (06-10, Wed, alone in its calendar week) -> 4, d 34 / d 35
(06-15 Mon / 06-16 Tue) -> median 34.5 -> 5.

Errors: mice in `data` missing from `meta` (lists them); `ref_col` or
`end_date` missing for a mouse with data -> `NA` week plus a warning naming
the mice; `data` already has any of the three output columns -> error
(name collision, like `meta_annotate()`).

### Cohort source (applies to 2 and 3) -- builder-check resolution

The cohort is the **meta workbook's `meta` sheet** (`meta_read(path)$meta`,
names cleaned by `janitor::clean_names()`: `mouse_ID`, `id`, `group`, `sex`,
`dob`, `bmt_date`, `death_date`, `early`, `reason`, ...), passed as
`workbook`. Reason: it holds all 23 enrolled mice; `dat$meta` lacks the 6
mice without `mouse_ID` *and* mice excluded mouse-wide (26-1-16 is dropped by
`db_load_experiment()`'s default). 26-1: all 6 non-DB mice have
`early = "yes"`, so "not in DB" is detected as `mouse_ID` NA or not in
`dat$meta` / exclusions, independent of `early`.

Future: Janis plans (separate session) to import the whole meta sheet into
the DB. So the cohort input is normalised internally: either workbook names
(`mouse_ID`, `id`, `death_date`, `dob`, `early`, `reason`) or DB names
(`mouse_id`, `end_date`, ...) are accepted, mapped to one internal shape. Then
the DB cohort can replace the workbook without an API change. *Claude's call.*

`workbook = NULL` falls back to `dat$meta` (+ exclusions' mouse ids) with a
message that mice without a DB record are not shown.

**Excluded** = excluded mouse-wide in every import unit, the same rule as
`fully_excluded_mice_()` (reuse it). Partial exclusions (one tissue / panel /
condition / file / date) do not make a mouse "excluded"; the grid shows them as
a note (`ok; excl. kidney (reason)`).

### 2. `meta_summary_table(workbook, dat = NULL, group_col = "group")`

File `R/meta_tables.R`. Returns a `gt` table: rows = characteristics, columns =
one per group level (factor order if `group` is a factor, else sorted), plus a
"not assigned" column only when an enrolled mouse has no group (26-1: the 2
removed at BMT).

Rows (counts are n, continuous are median (min-max)):
- Enrolled, n
- Sex, f / m
- Age at BMT, weeks (from `dob`; skipped if all NA)
- Strain / generation (distinct values joined; `line`/`generation`)
- Cages, n distinct
- BMT date (distinct dates)
- Dropped before / at BMT, n (`early` and (`bmt_date` NA or
  `death_date <= bmt_date`)) -- reasons as footnote
- Ended early after BMT, n (`early`, `death_date > bmt_date`) -- reasons
  footnote
- Excluded from analysis, n (+ mouse ids; mouse-wide rule above; from
  `dat$exclusions`, row omitted without `dat`)
- End-point, days post-BMT (mice not early)
- Analysed, n = mice with end-point data in at least one import unit and not
  excluded (mouse-wide). Early-ending mice with end-point data count
  (26-1-17 -> 26-1 analysed = 16). *Decided by Janis.* Needs the import
  units, so the function also takes optional `dat` (`meta_summary_table(workbook,
  dat = NULL)`, exclusions read from `dat$exclusions`); without `dat` the
  Analysed row is omitted.

### 3. `meta_mouse_table(dat, workbook = NULL)`

File `R/meta_tables.R`. Takes the whole `db_load_experiment()` list plus the
workbook sheet. One row per enrolled mouse (workbook rows; non-DB mice listed
by ear tag `id` with blank `mouse_id`), sorted by group then mouse id
(natural sort: 26-1-2 before 26-1-10). Returns a `gt` table.

Columns: `mouse_id`, `group`, `sex`, `cage`, age at BMT (weeks), end
(days post-BMT), then one column per import unit (names of `dat` except
`meta`, `exclusions` and elements that are not data frames with `mouse_id`).
Cell content:
- tissues present, when the unit has more than one distinct tissue
  (first letters, e.g. `K S L`);
- otherwise, when the unit has serial samples (any `collected_at` not equal to
  the mouse's `end_date`): `present/expected` where expected = distinct serial
  dates of the unit + 1 end-point (26-1 titer: complete = `4/4`, also for the
  early-ending 26-1-17 when its bleeds and end-point exist);
- otherwise a check mark.
Blank = no data. The `end_date` comes from `dat$meta` (or workbook
`death_date`). Last column `status`: `ok`, `early: <reason> (d <n>)`,
`excluded: <reason>`, `not in DB: <reason>`, plus partial-exclusion notes.
Rows not `ok` are highlighted.

### 4. `plot_timecourse(data, x = week_num, y = value, fill = group, subject = mouse_id, facet = NULL, scale = c("linear", "log"), colors = getOption("JanisHelpers.group_colors"), p_val = FALSE)`

File `R/plot_timecourse.R`. Input: long ELISA tibble after `meta_label_weeks()`
and `meta_annotate()` (has `group`).

- Rows with `NA` `y` are dropped with a message stating how many (no
  substitution; OOR< stays off the plot).
- Replicates are averaged first: mean `y` per `subject` x `x` x `fill`
  (x `facet`).
- Layers: thin per-mouse lines + small points (colour = `fill`, alpha ~0.4),
  thick group-median line + larger median points. No IQR band here (the
  individual lines show the spread; a band would interpolate between weeks --
  same reasoning as the summary plot). *Claude's call, changed from the first
  sketch.*
- X axis: when `x` is numeric, real time; breaks at the weeks present,
  labelled `"week <n>"`. When `x` is a factor (e.g. `week`), categorical.
- `lloq` column present and non-NA -> dashed horizontal line per facet (unique
  value; if several, the max). Absent / all-NA -> no line.
- `scale = "log"` uses `add_y_scale_()`; values must be > 0 (no
  pseudocount; error message if not).
- `facet` (bare, optional): `facet_wrap()` on it with free y.
- `p_val = TRUE`: per `x` (and facet) group comparison -- 2 groups Wilcoxon
  rank-sum (`rstatix`), > 2 Kruskal-Wallis; BH-adjusted across weeks within a
  facet; stars (`ns` omitted) above the week. Weeks where any group has < 2
  mice are skipped.
- Colours: the same `colors` / option logic as `plot_box()`, applied to both
  `colour` and `fill` (extend `add_fill_scale_()` or add a sibling for
  colour).

### 5. `plot_timecourse_summary(...)` -- same arguments as `plot_timecourse()` plus `min_n = 2`

Same file and preprocessing (shared internal for NA drop, replicate averaging,
axis, lloq, scale, stars, colours).

- Per group x week: median point, IQR error bar (25th-75th percentile,
  `stats::quantile(type = 7)`), median line joining weeks; groups dodged
  (`position_dodge(width = ~0.4 week)`).
- Points with n < `min_n` (26-1: week 4, irrel, n = 1) are drawn hollow,
  without error bar, and left out of the median line, which **bridges** them
  (irrel: week 3 -> week 5 directly). *Decided by Janis.*
- n per group printed inside the bottom of the panel under each week,
  coloured by group.
- Caption set by default: `"median and IQR; n per group below"` (user can
  override with `labs()`).

### 6. Wiring

- `analysis_26-1_db.Rmd` (OneDrive, outside the repo):
  - setup: `dat$workbook <- meta_read(...)$meta` (replaces `dat$dropouts`);
  - top "Cohort" section: `meta_summary_table(dat$workbook, dat)`;
  - QC: `meta_mouse_table(dat, dat$workbook)` replaces the
    `qc-overview` count, `qc-dropouts`, `qc-missing-samples`,
    `qc-exclusions` chunks;
  - titer chunk: `meta_label_weeks()` in the setup on `dat$titer`, then
    `plot_timecourse()` and `plot_timecourse_summary()`.
- `analysis-db` template: Overview section uses the two tables (dropouts
  optional, commented workbook read); ELISA section gets a commented
  timecourse example.
- CLAUDE.md data-structures map: one entry for the meta tables + weeks, one for
  the timecourse plots.

## Testing

- `test-meta_weeks.R`: the 26-1 date pattern above (synthetic tibbles) gives
  -1/1/3/4/5; pooling across a rounding boundary (d 31 / d 32 same calendar
  week -> one week); serial samples not pooled; missing-mouse error; collision
  error; ordered factor levels.
- `test-meta_tables.R`: counts in the summary with/without dropouts and
  exclusions (unmatched dropout counted once, matched one not double-counted);
  mouse grid presence cells and statuses; both return `gt_tbl`.
- `test-plot_timecourse.R`: returns `ggplot`; NA rows dropped with message;
  replicate averaging (layer data); lloq line present only with lloq; log with
  non-positive values errors; summary `min_n` hollow point / no errorbar;
  `p_val` builds without error on 2 groups.

## Left to the builder

- Exact `gt` styling, footnote wording, abbreviations.
- Internal helper names and how the colour scale helper is generalised.
- Whether star placement uses `ggpubr` or a hand-built `geom_text` (rstatix
  results).
