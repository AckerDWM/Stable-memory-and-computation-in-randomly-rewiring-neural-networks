# Generate a pool of synaptic weights sampled from the empirical distribution

set.seed(9812341)

P = function(s) {
  A = 100.7
  B = 0.02
  sigma1 = 0.022
  sigma2 = 0.018
  sigma3 = 0.150
  A*(1-exp(-s/sigma1))*(exp(-s/sigma2)+B*exp(-s/sigma3))
}

W = function(s) {
  (s/.2)*(s/(s+0.0314))
}

s = runif(1000000, 0, .2)
p = runif(1000000, 0, 23)
synaptic_size_pool = s[p<=P(s)]
synaptic_strength_pool = W(synaptic_size_pool)


# Load the pre-computed matrix of grid cell firing rates by running track position

load_input = function() {
  wd = getwd()
  setwd("/Users/danielacker/Google Drive/Memory-in-structurally-unstable-neural-networks/simulation_code/grid_cells/")
  grid = np$load("grid_cells-2d.npy") %>% t()
  setwd(wd)
  return(grid)
}