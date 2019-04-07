library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
np = reticulate::import("numpy")

source("simulation-code/models-1-and-2/synaptic-strength-pool.R")
source("simulation-code/models-1-and-2/model-2-methods.R")

grid = np$load("simulation-code/grid-cell-matrices/numpy-binaries/grid_cells-2d.npy") %>% t()

# simulate place cell activity in a rewiring network with competitive inhibition
sim = function(X, learning_rate, dropout_fun="E%-max") {
  # values to save
  place_cell_list = list()
  # parameters
  n_units = 2000 # number of place cells
  n_inputs = ncol(X) # number of grid cells
  # initialize place cells
  place_cells = initialize_place_cells(n_units, n_inputs)
  # day zero
  place_cells = calculate_firing_rates(place_cells=place_cells, X=X, dropout_fun=dropout_fun) # early phase
  place_cells = update_place_cell_weights(place_cells, X, learning_rate)
  place_cells = calculate_firing_rates(place_cells=place_cells, X=X, dropout_fun=dropout_fun) # late phase
  place_cell_list[[1]] = place_cells # save place cell data
  # days 1 to 60
  for (day in 1:60) {
    place_cells = turnover(place_cells) # synapse turnover
    place_cells = calculate_firing_rates(place_cells=place_cells, X=X, dropout_fun=dropout_fun) # early phase
    place_cells = update_place_cell_weights(place_cells, X, learning_rate)
    place_cells = calculate_firing_rates(place_cells=place_cells, X=X, dropout_fun=dropout_fun) # late phase
    place_cell_list[[day+1]] = place_cells # save place cell data
  }
  return(place_cell_list)
}

# set plasticity rates for simulations
eta = c(0, 1e-4) # extend this vector to add replicates

# run simulation and collect place field properties as a dataframe
## E%-max is the dropout function used in the paper
df_results = 
  sapply(eta, function(eta) {
    sim(grid, learning_rate=eta, dropout_fun="E%-max") %>%
      get_mean_place_field_properties() %>%
      mutate(Eta = eta)
  }, simplify=F) %>%
  bind_rows()

# generate a summary plot
ggplot(df_results, aes(Day, Median_drift, color=as.character(Eta))) +
  geom_point()
