readin_params <- function ( parfile = "inputs/par.csv" , cont_mat="inputs/cont_mat_2grp.csv" ){
# xx <- read.csv( "inputs/par.csv" )
xx = read.csv( parfile )
cont_mat = read.csv( cont_mat )
cont_vec = c( t( cont_mat ) )
params = list()
params = within( params, {
  ## set rng state
  seed <- xx[xx$par=="seed",]$val
  date_init <- as.Date("2020-02-07")
  tau <- xx[xx$par=="tau",]$val
  ndays <- as.integer(as.Date("2020-02-18") - as.Date("2020-02-07"))  + nrow(DAT)
  day_filter <- seq( 1, by=round(1/tau), length.out=ndays )
  nsteps <- round(ndays/tau) ## total number of steps
  delta1 <- xx[xx$par=="delta1",]$val# rate from infection to infectiousness
  delta2 <- xx[xx$par=="delta2",]$val# rate from infectiousness to symptom
  gamma <- xx[xx$par=="gamma",]$val # natural recovery rate
  frac_asymp <- xx[xx$par=="frac_asymp",]$val# fraction of infections that lead to asymptomatic infection
  time_intervention_start <- xx[xx$par=="time_intervention_start",]$val # as.Date("2020/02/22") # highest level
  time_intervention_stop <- xx[xx$par=="time_intervention_stop",]$val # as.Date("2020/02/22") # highest level
  cm <- cont_vec
  rate_isol <- xx[xx$par=="rate_isol",]$val # delay until isolation for high risk
  rel_rate_isol_asymp <- xx[xx$par=="rel_rate_isol_asymp",]$val # delay until isolation for asymptomatically infected
  max_rate_isol <- xx[xx$par=="max_rate_isol",]$val # delay before isolation for high risk
  R01 <- xx[xx$par=="R01",]$val
  R02 <- xx[xx$par=="R02",]$val
  final_R01 <- xx[xx$par=="final_R01",]$val
  final_R02 <- xx[xx$par=="final_R02",]$val
  #
  ## initial conditions, list within list
  ## use within() to modify empty list, as above
  pop_high <- round( xx[xx$par=="pop_high",]$val )
  pop_low <- round( xx[xx$par=="pop_low",]$val )
  ngrroups <- xx[xx$par=="ngroups",]$val
  
  pop <- c( pop_high, pop_low )
  
  I01 <- xx[xx$par=="I01",]$val
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
return( list(params=params, savepar=xx, cont_mat=cont_mat) )
}