library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)
np = reticulate::import("numpy")

source("simulation-code/models-1-and-2/synaptic-strength-pool.R")

grid = np$load("simulation-code/grid-cell-matrices/numpy-binaries/grid_cells-2d.npy") %>% t()

simulate = function(
  X, n_inputs, n_outputs, turnover, eta, pretraining=T, post_training=T) {
  
  threshold_function = function(x) {x[x - quantile(x, .9) < 0]=0;x}
  
  # Initialize inputs and synaptic weights
  orginal_input_idx = sample(1:10000, n_inputs)
  new_input_idx = sample(which(!(1:10000) %in% orginal_input_idx), turnover)
  input = X[,orginal_input_idx]
  weights = replicate(n_outputs, sample(synaptic_strength_pool, size=n_inputs, replace=T))
  replacement_input = X[,new_input_idx]
  
  # save original weights
  df_weights = data.frame(original_weights = c(weights))
  
  # Activate the output units pre-turnover
  response = input %*% weights
  response %<>% apply(2, threshold_function)
  
  # Pre-training
  if (pretraining) {
    dw = np$dot(t(input), response)
    weights = weights + eta*dw
    weights %<>% apply(2, function(x) 148.7*x/sum(x))
  }
  
  # Save trained weights
  df_weights$trained_weights = c(weights)
  
  # Save weights and inputs
  inputs_pre_turnover = input
  weights_pre_turnover = weights
  
  # Synpase turnover
  input[,1:turnover] = replacement_input
  weights[1:turnover,] = sample(synaptic_strength_pool, size=turnover*n_outputs, replace=T)
  
  # Activate the output units post-turnover
  response = input %*% weights
  response %<>% apply(2, threshold_function)
  
  # Post-training
  if (post_training) {
    dw = np$dot(t(input), response)
    weights = weights + eta*dw
    weights %<>% apply(2, function(x) 148.7*x/sum(x))
  }
  
  # Save trained/post-turnover weights
  df_weights$post_turnover_weights = c(weights)
  
  # Sum synaptic potentials across replaced inputs
  PSP_pre_turnover = inputs_pre_turnover[,1:turnover] %*% weights_pre_turnover[1:turnover,]
  PSP_post_turnover = input[,1:turnover] %*% weights[1:turnover,]
  
  # Get correlation between PSPs pre and post-turnover
  correlation = cor(c(PSP_pre_turnover), c(PSP_post_turnover))
  
  list(w=df_weights, c=correlation)
}

# perform simulations with increasing rates of turnover
# expressed as the number of synapse replaced in one day
turnover = 1200*seq(.1, 1, .1)
result_varying_turnover = 
  sapply(turnover, function(turnover) {
    replicate(100, {
      simulate(
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
      simulate(
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
