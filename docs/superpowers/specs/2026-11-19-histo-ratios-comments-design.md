# Histology prep: multiple ratios, comment column, slimmer QC report

Short design (small change). Builds on
`2026-11-18-elisa-histo-prep-qc-export-design.md`.

## Intent

New columns in the curated histology sheet (e.g. `abnormal_gloms`) already
flow through read -> QC -> export as metrics. The user wants (1) several
percent ratios in one step (`pct_cresc`, `pct_necrosis`, `pct_abnormal`),
(2) a free-text comment column that does not break the read, and (3) a
shorter QC report.

## Decisions

1. **Vectorised `histo_calc_ratio()`.** `numerator` may be a named character
   vector; names are the new metric names (`name` is then not used).
   `denominator` is length 1 (recycled) or one per numerator. The old
   scalar call `histo_calc_ratio(data, num, den, name)` keeps working.
   Same validation as today (new names must not exist / be unique, every
   referenced metric must exist). The prep template uses the vector form.
2. **Comment column (histo only).** A `values` column named `comment` or
   `comments` (case-insensitive) is read as text, not a metric. It is
   carried as a `comment` column on the long table (repeated per metric
   row, `NA` when empty). Not exported to the database -- the curated
   `.xlsx` stays the record; `db_export_histo()` ignores the column. The QC
   report gets a "Comments" section: one row per mouse x tissue x date with
   a non-empty comment, flagged when excluded. The blank template and its
   README mention the column.
3. **QC report.** Remove the "By batch" and "Over time" outlier sections.
   `batch_col` stays a prep setting (overview row, shared with FACS/ELISA).

## Approaches considered

- Ratios: separate `histo_calc_ratios()` vs vectorising the existing
  function -- chose vectorising (one function to know, old calls intact).
- Comments: report-only (chosen) vs a new `histo_comments` DB table --
  deferred; no schema change now.

## Non-goals

ELISA comment handling; auto-skipping ratios / rules whose metrics are absent.
