# article_bigsimr
Rmarkdown article related to the bigsimr package

In order to build this article from scratch on your own PC, please install the dependencies by running the following:
  
  ```r
renv::restore()
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

## Generate Data and Figures

```r
source("setup.R") # this will take a while
```
