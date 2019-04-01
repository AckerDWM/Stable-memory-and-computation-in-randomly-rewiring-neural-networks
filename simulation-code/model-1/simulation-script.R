source("requirements.R")
source("methods.R")
source("sim-EPSC-cor.R")

all_grid_cells = load_input()


# simulate stability at varying turnover rates

set.seed(9451243)

turnover = 1200*seq(.1, 1, .1)

result_varying_turnover = sapply(turnover, function(turnover) {
  replicate(
    100, 
    {
      simulate(
        n_cosines=1200, n_outputs=1, turnover=turnover, 
        eta=1e-4, pretraining=T, post_training=T
      )$c
    }
  )
})

colnames(result_varying_turnover) = turnover

result_varying_turnover %<>%
  melt()%>% 
  rename(Replicate = Var1, `Turnover` = Var2, `EPSC correlation` = value)


# simulate omission of session one learning

set.seed(9532122)

pretraining = c(T, F)

result_no_session_one = sapply(pretraining, function(pretraining) {
  replicate(
    100, 
    {
      simulate(
        n_cosines=1200, n_outputs=1, turnover=120, 
        eta=1e-4, pretraining=pretraining, post_training=T
      )$c
    }
  )
})

colnames(result_no_session_one) = pretraining
result_no_session_one %<>%
  melt() %>%
  mutate(Var2 = factor(Var2, levels=c("TRUE", "FALSE"))) %>%
  mutate(Var2 = plyr::revalue(Var2, replace=c(
    "TRUE"="LTP", "FALSE"="No LTP"
  ))) %>% 
  rename(Replicate = Var1, Condition = Var2, `EPSC correlation` = value)