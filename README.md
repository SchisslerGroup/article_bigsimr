# article_bigsimr
Rmarkdown article related to the bigsimr package

In order to build this article from scratch on your own PC, please install the dependencies by running the following:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

devtools::install_deps()
devtools::install_version("TCGA2STAT", version = "1.2")
devtools::install_github("SchisslerGroup/bigsimr", ref = "dev-julia")
```

## Other Setup

During setup, `JuliaCall` will try to look for your julia installation by first checking the environment variable `JULIA_HOME`. If it is not set, then it will try to automatically install julia. This can sometimes be problematic. To get around this, first install julia, then edit your `.Renviron` file to include the following:

```
JULIA_HOME = "/path/to/julia/bin/"
```

You can edit your `.Renviron` file by calling

```r
usethis::edit_r_environ("user")
```
