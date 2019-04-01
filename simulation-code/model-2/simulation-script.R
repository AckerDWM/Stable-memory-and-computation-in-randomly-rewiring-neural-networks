source("requirements.R")
source("methods.R")
source("sim.R")

set.seed(83721)

# set learning rates for simulations
eta = c(0, 1e-4)

# run simulation and collect place field properties as a dataframe
df_results = 
  sapply(eta, function(eta) {
    sim(learning_rate=eta, dropout_fun="E%-max") %>%
      get_mean_place_field_properties() %>%
      mutate(Eta = eta)
  }, simplify=F) %>%
  bind_rows()