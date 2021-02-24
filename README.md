# article_bigsimr
Rmarkdown article related to the bigsimr package

In order to build this article from scratch on your own PC, please install the dependencies by running the following:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

devtools::install_deps()
devtools::install_version("TCGA2STAT", version = "1.2")
devtools::install_github("SchisslerGroup/bigsimr", ref = "dev-julia")
```
