devtools::load_all()
library(tidyverse)

dat_long <- benchmark_dependences %>%
  select(-total_time, -n_sim, -needed_near_pd) %>%
  pivot_longer(cols = c(corr_time:sim_time), names_to = "Step", values_to = "Time") %>%
  mutate(dim = factor(dim),
         Step = factor(Step,
                       levels = c("corr_time",
                                  "adj_time",
                                  "admiss_time",
                                  "sim_time"),
                       labels = c("Compute Correlation",
                                  "Adjust Correlation",
                                  "Check Admissibility",
                                  "Simulate Data"))) %>%
  rename(Correlation = corr_type, Dimensions = dim)

dat_long %>%
  filter(Dimensions %in% c(100, 250, 500, 1000)) %>%
  ggplot(aes(Correlation, Time, fill=Step)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Time (Seconds)") +
  facet_grid(. ~ Dimensions, scales = "free") +
  theme(axis.text.x = element_text(angle = -90))


dat_long %>%
  filter(Dimensions %in% c(2500, 5000, 10000)) %>%
  mutate(`Time (minutes)` = Time / 60) %>%
  ggplot(aes(Correlation, `Time (minutes)`, fill=Step)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 30, 45),
                     limits = c(0, 30),
                     minor_breaks = NULL) +
  facet_grid(. ~ Dimensions) +
  theme(axis.text.x = element_text(angle = -90))


