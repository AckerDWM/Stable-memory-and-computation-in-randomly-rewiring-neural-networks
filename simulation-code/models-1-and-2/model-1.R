library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)
np = reticulate::import("numpy")

source("simulation-code/models-1-and-2/synaptic-strength-pool.R")
source("simulation-code/models-1-and-2/model-1-simulate_epsc.R")

grid = np$load("simulation-code/grid-cell-matrices/numpy-binaries/grid_cells-2d.npy") %>% t()

# perform simulations with increasing rates of turnover
# expressed as the number of synapse replaced in one day
turnover = 1200*seq(.1, 1, .1)
result_varying_turnover = 
  sapply(turnover, function(turnover) {
    replicate(100, {
      simulate_epsc(
        X=grid, n_inputs=1200, n_outputs=1, turnover=turnover, 
        eta=1e-4, pretraining=T, post_training=T
      )$c
    })
  }) %>%
  { colnames(.) = turnover;. } %>% 
  melt() %>% 
  rename(Replicate=Var1, `Turnover`=Var2, `EPSC correlation`=value)

# perform simulations omiting learning on first session
# mimicking NMDAR blockade
pretraining = c(T, F)
result_no_session_one = 
  sapply(pretraining, function(pretraining) {
    replicate(100, {
      simulate_epsc(
        X=grid, n_inputs=1200, n_outputs=1, turnover=120, 
        eta=1e-4, pretraining=pretraining, post_training=T
      )$c
    })
  }) %>%
  { colnames(.)=pretraining;. } %>%
  melt() %>%
  mutate(Var2 = factor(Var2, levels=c("TRUE", "FALSE"))) %>%
  mutate(Var2 = plyr::revalue(Var2, replace=c(
    "TRUE"="LTP", "FALSE"="No LTP"
  ))) %>% 
  rename(Replicate = Var1, Condition = Var2, `EPSC correlation` = value)

# plot results
g_varying_turnover = 
  result_varying_turnover %>%
  ggplot() +
  aes(factor(Turnover), `EPSC correlation`) +
  geom_boxplot()

g_no_sessin_one = 
  result_no_session_one %>%
  ggplot() +
  aes(Condition, `EPSC correlation`) +
  geom_boxplot()

plot_grid(g_varying_turnover, g_no_sessin_one)
