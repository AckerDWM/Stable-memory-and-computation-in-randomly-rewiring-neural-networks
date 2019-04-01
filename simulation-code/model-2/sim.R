# simulate place cell activity in a rewiring network with competitive inhibition
sim = function(learning_rate, dropout_fun="E%-max", grid_type="random") {
  # values to save
  place_cell_list = list()
  # parameters
  X = load_input(grid_type) # grid cell firing rates
  n_units = 2000 # number of place cells
  n_inputs = 10000 # number of grid cells
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
