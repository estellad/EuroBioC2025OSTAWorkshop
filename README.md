# Workflow: Visium HD cell-level

Authors: Yixing E. Dong[^readme-1], Ellis Patrick[^readme-2].

[^readme-1]: University of Lausanne, Lausanne, Switzerland

[^readme-2]: University of Sydney, Sydney, Australia

## Instructor name and contact information

-   Yixing Estella Dong ([estella.yixing.dong\@gmail.com](mailto:estella.yixing.dong@gmail.com){.email})
-   Ellis Patrick ([ellis.patrick\@sydney.edu.au](mailto:ellis.patrick@sydney.edu.au){.email})

## Workshop Description

In this instructor-led live demo, we analyse Visium HD data segmented to cells, demonstrating use of `SpatialExperiment` and `sf` classes in R to import, organise, quality control, clustering, marker gene identificaiton, and label transfer to annotate the data. Spatial transcriptomics data with cell types annotated can reveal key insights with neighbourhood enrichment and spatial statistics, which will be demonstrated in a second half of the workshop. The complete analysis offered by existing tools highlights the ease with which researchers can turn the raw counts from a Visium HD experiment into biological insights using Bioconductor and CRAN. The complete workflow is available at <https://lmweber.org/OSTA/>

## To run locally.

Clone the repo and follow the vignette in the `vignettes` folder. You will also need to install all the packages needed for the workshop. This can be done by

`devtools::install_github("estellad/EuroBioC2025_OSTA_Workshop")`
