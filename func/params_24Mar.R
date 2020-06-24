params <- list()
params <- within( params, {
  ## set rng state
  seed <- 0
  date_init <- as.Date("2020/02/07")
  tau <- 0.01 # in days
  ndays <- as.integer(as.Date("2020/02/19") - as.Date("2020/02/07") + 1) + length(DAT$case_total)
  day_filter <- seq( 1, by=round(1/tau), length.out=ndays )
  nsteps <- round(ndays/tau) ## total number of steps
  delta <- 1/5.2 # 1/delta=mean incubation period
  gamma <- 0 # natural recovery rate
  frac_asymp <- 0.0 # fraction of infections that lead to asymptomatic infection
  time_intervention_start <- 16 # as.Date("2020/02/22") # highest level 
  cm <- c( 0.999, 0.001, 0.001, 0.999 )
  rate_isol_x1 <- 1/4 # delay until isolation for high risk
  rate_isol_x2 <- 1/4 # delay until isolation for high risk
  rel_rate_isol_x1_asymp <- 1.0 # delay until isolation for asymptomatically infected
  rel_rate_isol_x2_asymp <- 1.0 # delay until isolation
  daily_frac_reduc1 <- 0.1 # prop reduction of gamma(i.e, isolation rate increase) 
  daily_frac_reduc2 <- 0.1 # prop reduction of gamma(i.e, isolation rate increase) 
  min_delay1 <- 0.5 # delay before isolation for high risk
  min_delay2 <- 0.5 # delay before isolation for low risk  
  beta1 <- 2
  beta2 <- 0.15
  min_beta1 <- 0.2 #
  min_beta2 <- 0.2 # 
  ## initial conditions, list within list
  ## use within() to modify empty list, as above
  # pop_high <- POP_SHINCHEONJI_DAEGU*AVG_HOUSEHOLD_SIZE_DAEGU
  pop_high <- POP_SHINCHEONJI_DAEGU
  pop_low <- POP_DAEGU - pop_high 
  pop <- c( pop_high, pop_low )
  I01 <- 10
  init <- within(list(), {
    S1 <- pop[1]-I01
    E1_1 <- 0
    E1_2 <- 0
    I1 <- I01
    A1 <- 0
    R1 <- 0
    X1 <- 0
    CI1 <- 0
    S2 <- pop[2]
    E2_1 <- 0
    E2_2 <- 0
    I2 <- 0
    A2 <- 0
    R2 <- 0
    X2 <- 0
    CI2 <- 0
  })
})
