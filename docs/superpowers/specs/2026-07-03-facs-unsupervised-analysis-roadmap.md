# Roadmap: unsupervised analysis of FACS data

**Date:** 2026-07-03
**Status:** Big-picture plan — not yet a stage-level design. Revisit this
file at the start of each future session on this initiative before drafting
a stage's own design doc.

---

## Goal

Add unsupervised, cell-level analysis of flow cytometry data to the
package: cluster raw single-cell events (not just FlowJo's manually gated
populations) to discover phenotypes, visualize them, and test whether
cluster abundance differs between experimental groups.

This is a substantially larger scope than the existing `facs_` functions
(which all operate on FlowJo's already-gated summary statistics). It is
being decomposed into four sequential stages, each with its own
spec -> plan -> implementation cycle. Stages are sequential, not
independent — each reads the output of the one before it.

---

## Decisions made so far (brainstorming session, 2026-07-03)

- **Analysis level:** cell-level (raw single-cell events), not sample-level
  summary-statistic clustering.
- **Scope:** full pipeline — cluster, visualize, annotate, test differential
  abundance.
- **Raw data availability:** both `.fcs` and `.wsp` files are available
  side by side.
- **Pre-filtering:** replay the existing FlowJo `.wsp` gating tree against
  the raw `.fcs` events (via `CytoML`/`flowWorkspace`, both new Bioconductor
  dependencies) rather than a lightweight manual threshold filter or no
  filter at all. Reuses the existing gating investment and keeps results
  comparable to already-gated data.
- **Starting gate for clustering:** not fixed to a convention (e.g. always
  "live CD45+ singlets") — must be a caller-supplied parameter, since it
  will vary by experiment/panel.
- **Marker/channel selection for clustering:** explicit, caller-supplied
  list per call — no auto-detection of "all markers except gating
  markers."
- **Clustering algorithm:** FlowSOM (self-organizing map + consensus
  metaclustering). Chosen over PhenoGraph for speed at scale and because
  cluster granularity can be picked after the fact via metaclustering.
- **Dimensionality reduction:** UMAP (over tSNE) — faster, better global
  structure preservation, current default in modern cytometry pipelines.
- **Differential abundance testing:** diffcyt-style GLMM (Bioconductor
  `diffcyt`, likely pulling in `edgeR`/`lme4`), not a reuse of the existing
  `rstatix`-based `analysis_` domain. Deliberate trade-off: this is a
  second, separate statistical framework alongside the package's existing
  `gtsummary`/`rstatix` one, but purpose-built for count-based cluster
  abundance data (overdispersion, repeated measures/batch effects) in a way
  the existing Wilcoxon/Kruskal machinery is not.

---

## Stage breakdown

### Stage 1 — Raw FCS import + gate replay (not yet designed)
Read raw single-cell events from `.fcs` files, filtered down to an
arbitrary gate path defined in the matching `.wsp` workspace (via
`CytoML`/`flowWorkspace`), restricted to a caller-supplied list of
markers/channels. Foundation for every later stage — brainstorm this one
next.

### Stage 2 — FlowSOM clustering (not yet designed)
Run FlowSOM (SOM + consensus metaclustering) on Stage 1's output. Produces
per-event cluster assignments and per-sample cluster frequencies.

### Stage 3 — UMAP dimensionality reduction (not yet designed)
Compute UMAP embeddings on the gated events (or a sample of them) for
visualization, likely colored/faceted by Stage 2's cluster assignments.

### Stage 4 — Differential abundance testing (not yet designed)
Test whether Stage 2's per-sample cluster frequencies differ between
experimental groups, using a diffcyt-style GLMM. Will need to integrate
with the existing `meta_` domain for group/design information.

---

## Optional additions

### Audit trail (optional, not yet designed)
Each wrapper function in this pipeline (`facs_read_fcs_gated()`,
`facs_cluster_flowsom()`, `facs_reduce_umap()`, and later the Stage 4
differential-abundance function) takes parameters that materially affect
the result (gate path, markers, grid/metacluster counts, seeds, UMAP
neighbors/min_dist, etc.). To keep results reproducible, add an opt-in
audit trail: a dated document (one per run/session) recording every call
to these wrapper functions with its full argument list, so a later reader
can reconstruct exactly how a given result was produced without
re-deriving it from script history. Open questions to settle if/when this
is brainstormed: where the log lives (per-project file vs. per-call
sidecar), what triggers a new dated document vs. appending to an existing
one, and whether this is a cross-cutting helper (its own domain, e.g.
`wrangle_` or a new `audit_` prefix) or logic embedded in each wrapper.

---

## Open questions for future stages

- Which domain prefix these new functions should use (`facs_` fits the
  cytometry-specific stages; the diffcyt-based stage may sit better in
  `analysis_`, or may need to stay in `facs_` for consistency with the rest
  of the pipeline — decide when Stage 4 is brainstormed).
- New dependencies required (`CytoML`, `flowWorkspace`, `FlowSOM`,
  `uwot` or similar for UMAP, `diffcyt` and its own dependencies) are all
  Bioconductor/CRAN packages not yet in `DESCRIPTION` — each stage's design
  should confirm exact package names and add them to `Imports` at
  implementation time.
