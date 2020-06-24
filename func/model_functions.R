run_model_tauleap <- function( params ){
  
  out <- seir_2grp_tauleap( params )
  day_filter <- seq( 1, by=round(1/params$tau), length.out=params$ndays )
  out <- out[ day_filter, ]
  
  n2 <- nrow(out)
  n1 <- length(DAT$daily_total) - 1
  
  num_case_high <- out$X1[(n2-n1):n2] - out$X1[(n2-n1-1):(n2-1)]
  num_case_low <- out$X2[(n2-n1):n2] - out$X2[(n2-n1-1):(n2-1)]

  return( list( num_case_high, num_case_low ) )
}


calc_distance <- function( Dstar, Dstar2 ){
  
  id <- 12:18  # spiky area, sum over time will be used to measure the difference
  id2 <- 10:11
  # dist_inc_high <- sqrt( sum( (Dstar - DAT["daily_shincheonji"])^2, na.rm=T ) )
  # dist_inc_low <- sqrt( sum( (Dstar2 - DAT["daily_low"])^2 , na.rm=T ) )
  dist_inc_high <- sqrt( sum( (Dstar[-id] - DAT$daily_shincheonji[-id])^2, na.rm=T ) + (sum(Dstar[id])-sum(DAT$daily_shincheonji[id]))^2 )
  dist_inc_low <- sqrt( sum( (Dstar2[-c(id,id2)] - DAT$daily_low[-c(id,id2)])^2 , na.rm=T ) + (sum(Dstar2[id])-sum(DAT$daily_low[id]))^2 + (sum(Dstar2[id2])-sum(DAT$daily_low[id2]))^2 )
  # dist_inc_low <- sqrt( sum( (Dstar2[-id] - DAT$daily_low[-id])^2 , na.rm=T ) + (sum(Dstar2[id])-sum(DAT$daily_low[id]))^2 )
  
  
  return( c( dist_inc_high, dist_inc_low) )
}


# Perturbation kernel 
rK <- function( mean, sigma ){   
  return( rtmvnorm( 1, mean=mean, sigma=sigma, lower=lb, upper=ub) ) 
}

#  Identity function: H(x)= 1 if x=T
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior_non_zero <- function(par){
  prod( sapply(1:length(par), function(a) H(par[a]-lb[a]) * H(ub[a]-par[a]) ) )
}

