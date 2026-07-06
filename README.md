
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JanisHelpers

<!-- badges: start -->
<!-- badges: end -->

Einzelne kleine Funktionen für alltägliche Probleme

- report_knit_dated: put in yaml behind "knit:" to automatically knit
  dated versions to specified directory
- Auto-knit-report: rmarkdown template das automatisch datiert dateien
  in spezifische subdirectories erstellt
- protocol-template: template für Protokolle
- Import FowJO workspace
- Automatically import and save FlowJo workspace to Excel

The unsupervised FACS pipeline (`facs_read_fcs_gated()`, `facs_cluster_flowsom()`,
`facs_test_cluster_abundance()`) needs `CytoML`/`flowWorkspace`/`flowCore`,
`FlowSOM`, and `diffcyt`/`SummarizedExperiment`, which are Bioconductor
packages not available on CRAN. Install them first:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("CytoML", "flowWorkspace", "flowCore", "FlowSOM", "diffcyt", "SummarizedExperiment"))
```

Updaten des privaten Pakets:
"devtools::install_github("JanisSonnemann/JanisHelpers", auth_token =
gh::gh_token())"

Erstellung eigener Pakete über: -
"<https://r-pkgs.org/whole-game.html#create_package>" -
"<https://r-pkgs.org/code.html>"
