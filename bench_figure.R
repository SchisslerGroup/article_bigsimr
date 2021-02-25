devtools::load_all()
library(tidyverse)

dat_long <- benchmark_dependences %>%
  select(-total_time, -n_sim, -needed_near_pd) %>%
  pivot_longer(cols = c(corr_time:sim_time), names_to = "Step", values_to = "Time") %>%
  mutate(dim = factor(dim),
         Step = factor(Step,
                       levels = c("corr_time", "adj_time", "admiss_time", "sim_time"),
                       labels = c("Compute Correlation",
                                  "Adjust Correlation",
                                  "Check Admissibility",
                                  "Simulate Data"))) %>%
  rename(Correlation = corr_type, Dimensions = dim)

ggplot(dat_long, aes(Dimensions, Time, fill=Step)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Time (Seconds)") +
  facet_grid(Correlation ~ ., scales = "free")
