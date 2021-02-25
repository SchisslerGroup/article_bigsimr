source("renv/activate.R")

message("Setting up Julia")
# JuliaCall::install_julia()

if (!dir.exists("data"))
  dir.create("data")

if (!dir.exists("results"))
  dir.create("results")

if (!dir.exists("fig"))
  dir.create("fig")

message("Generating Data")
if (!file.exists("data/example_brca.rda")) {
  source("data-raw/brca.R", local = TRUE)
}



message("Generating Figures")
