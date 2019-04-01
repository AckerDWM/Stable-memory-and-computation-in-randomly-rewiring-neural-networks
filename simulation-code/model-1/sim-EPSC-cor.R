# simulation function

simulate = function(
  n_cosines, n_outputs, turnover, eta, pretraining=T, post_training=T) {
  
  threshold_function = function(x) {x[x - quantile(x, .9) < 0]=0;x}
  
  # Initialize inputs and synaptic weights
  orginal_input_idx = sample(1:10000, n_cosines)
  new_input_idx = sample(which(!(1:10000) %in% orginal_input_idx), turnover)
  input = all_grid_cells[,orginal_input_idx]
  weights = replicate(n_outputs, sample(synaptic_strength_pool, size=n_cosines, replace=T))
  replacement_input = all_grid_cells[,new_input_idx]
  
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
  
  return(list(w=df_weights, c=correlation))
}
