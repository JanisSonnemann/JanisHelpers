# QC reports: mark excluded samples on point plots

## Intent

When reviewing a prep QC report, see at a glance whether a flagged outlier is
already explained by the meta workbook's `exclusions` sheet.

## Goals / non-goals

- Goal: every plot with one point per sample marks excluded samples as a
  cross (`x`, shape 4); others stay filled circles. Legend entry "excluded"
  appears only when the plot has at least one excluded point.
- Non-goal: changing statistics. Excluded samples still feed the boxes,
  whiskers and outlier flags/labels exactly as now.
- Non-goal: count/bar plots (ELISA empty-value and range charts).
- Anti-goal: hiding or dropping excluded points.

## Scope (plots)

| Report | Plot | Helper | Excluded source |
|---|---|---|---|
| FACS | individuals (hierarchy, counts, MFI) | `qc_plot_individuals_()` | `scored$excluded` (file-based, already there) |
| FACS | batch plots | `qc_plot_batch_()` | `scored$excluded` |
| ELISA / histo | outlier plot (group and batch view) | `qc_table_plot_()` | `excluded_keys` via `qc_is_excluded_()` |
| ELISA | replicate-CV plot | `elisa_qc_cv_plot_()` | `excluded_keys` via `qc_is_excluded_()` |

## Approaches

1. **Cross shape (picked).** Shape aesthetic; works on batch plots where
   colour = batch; a labelled cross reads as "explained outlier".
2. Ring overlay -- clutters dense facets, awkward legend.
3. "(excl.)" in labels -- misses excluded points inside the whiskers.

## Design

- Shared helper in `R/facs_qc.R`: a `scale_shape_manual()` mapping
  `FALSE` -> 16 (filled circle), `TRUE` -> 4 (cross), legend title
  "excluded", only the `TRUE` key shown (`breaks = TRUE`).
- Plot helpers map `shape = .excluded` on `geom_point()`; jitter seed and
  label layers unchanged, so labels still sit on their points.
- `qc_table_plot_()` and `elisa_qc_cv_plot_()` gain `excluded_keys = NULL`;
  `NULL` / no matches = all circles (current look, older callers unaffected).
  FACS helpers read the existing `excluded` column (missing column = FALSE).
- ELISA / histo templates pass `excluded_keys` to these calls; report text
  above the plots mentions "crosses = on the exclusions sheet".
- Tests: excluded points get shape 4, others 16, in each helper.

## Left to the builder

Exact wording of the report notes; whether the cross is slightly larger than
the circles for visibility.

## Addendum 2026-11-20: label colours

On plots that label points (FACS individuals, ELISA / histology outlier, ELISA
CV), excluded samples get a red `mouse_ID` label on every point; outliers
(and high-CV samples) get a black one. An excluded outlier gets only the red
label. Two fixed-colour label layers (`qc_label_layers_()`), so no colour
scale is touched. Batch plots stay unlabelled.
