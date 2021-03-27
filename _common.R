options(digits = 4)

knitr::opts_chunk$set(
  # Text Results
  eval = TRUE,
  echo = FALSE,
  results = 'markup',
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  # Code Decoration
  comment = NA,
  background = '#F7F7F7',
  # Cache
  cache = TRUE,
  # Plots
  fig.show = "hold",
  fig.width = 6.5,
  fig.asp = 0.618,  # 1 / phi
  out.width = "85%",
  fig.align = 'center'
)

# options(dplyr.print_min = 6, dplyr.print_max = 6)

seed_number <- 2021-03-23
set.seed(seed = seed_number)
