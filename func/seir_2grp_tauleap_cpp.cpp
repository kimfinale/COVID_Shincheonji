#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


List seir_2grp_tauleap(List params) {
// chained operations are tricky in cpp
// pull out list w/in list into its own object
List init = params["init"];

// use Rcpp as() function to "cast" R vector to cpp scalar
int nsteps = as<int>(params["nsteps"]);

// ints (IntegerVector), since rpois returns double,

NumericVector SS1( nsteps, init["S1"] );

NumericVector EE1_1( nsteps, init["E1_1"] );
NumericVector EE1_2( nsteps, init["E1_2"] ); // infectious, but did not developed symptoms
NumericVector II1( nsteps, init["I1"] ); // infectious and symptomatic
NumericVector AA1( nsteps, init["A1"] ); //asymptomatic
NumericVector XX1( nsteps, init["X1"] ); //isolated
NumericVector RR1( nsteps, init["R1"] );
NumericVector CI1( nsteps, init["CI1"] ); // cumulative

NumericVector SS2( nsteps, init["S2"] );
NumericVector EE2_1( nsteps, init["E2_1"] );
NumericVector EE2_2( nsteps, init["E2_2"] );
NumericVector II2( nsteps, init["I2"] );
NumericVector AA2( nsteps, init["A2"] );
NumericVector XX2( nsteps, init["X2"] );
NumericVector RR2( nsteps, init["R2"] );
NumericVector CI2( nsteps, init["CI2"] ); // cumulative
// fill time w/zeros
NumericVector time(nsteps);

// pull out params for easy reading 
NumericVector cm = params["cm"];//contact matrix
double delta1 = params["delta1"];
double delta2 = params["delta2"];
double gamma = params["gamma"];
double tau = params["tau"];
double fa = params["frac_asymp"];
// double rrxa = params["rel_rate_isol_asymp"]; //relative isolation rate for asymptomatic patients
double time_intervention_start = params["time_intervention_start"];
double time_intervention_stop = params["time_intervention_stop"];
double dur_intervention = time_intervention_stop - time_intervention_start;
double max_rate_x = params[ "max_rate_isol" ];
double final_R01 = params[ "final_R01" ];
double final_R02 = params[ "final_R02" ];
double R01 = params[ "R01" ];
double R02 = params[ "R02" ];

// Rprintf("the value of beta[1] : %f \n", beta1);
// Rprintf("the value of beta[2] : %f \n", beta2);

// Calculate the number of events for each step, update state vectors
for( int istep = 0; istep < (nsteps-1); istep++ ) {
  // pull out this step's scalars for easier reading
  // and to avoid compiler headaches
  double iS1 = SS1[istep]; //  rpois returns double
  double iS2 = SS2[istep];
  double iE1_1 = EE1_1[istep];
  double iE1_2 = EE1_2[istep];
  double iE2_1 = EE2_1[istep];
  double iE2_2 = EE2_2[istep];
  double iI1 = II1[istep];
  double iI2 = II2[istep];
  double iA1 = AA1[istep];
  double iA2 = AA2[istep];
  double iX1 = XX1[istep];
  double iX2 = XX2[istep];
  double iR1 = RR1[istep];
  double iR2 = RR2[istep];
  double iCI1 = CI1[istep];
  double iCI2 = CI2[istep];
  
  /////////////////////////
  // State Equations
  /////////////////////////
  
  // R::rpois always returns a single value
  // to return multiple (e.g. Integer/NumericVector, 
  // use Rcpp::rpois(int ndraw, param) and friends
  // Prevent negative states
  double rate_x = params["rate_isol"]; //initialize
  double rrxa = params["rel_rate_isol_asymp"]; //relative isolation rate for asymptomatic patients
  double beta1 = R01 / (fa/(gamma+rate_x*rrxa) + (1-fa)/(gamma+rate_x));
  double beta2 = R02 / (fa/(gamma+rate_x*rrxa) + (1-fa)/(gamma+rate_x));
  double final_beta1 = final_R01 / (fa/(gamma+rate_x*rrxa) + (1-fa)/(gamma+rate_x));
  double final_beta2 = final_R02 / (fa/(gamma+rate_x*rrxa) + (1-fa)/(gamma+rate_x));

  double time_diff = (istep+1)*tau - time_intervention_start;
  
  if( 0 < time_diff && time_diff <= dur_intervention ){
    rate_x = rate_x + time_diff*(max_rate_x - rate_x)/dur_intervention;
    beta1 = beta1 - time_diff*(beta1 - final_beta1)/dur_intervention;
    beta2 = beta2 - time_diff*(beta2 - final_beta2)/dur_intervention;
    rrxa = time_diff/dur_intervention;
  }
  
  else if( dur_intervention < time_diff ){
    beta1 = final_beta1;
    beta2 = final_beta2;
    rrxa = 1;
    rate_x = max_rate_x;
  }
  
  if( (istep+1)*tau < (time_intervention_start-3) ){// case identification occurs from Feb 18
    rate_x = 0;
  };
  
  double I1 = iI1 + iA1 + iE1_2; // asymptomatic individuals have the same transmissibility
  double I2 = iI2 + iA2 + iE2_2;
  double N1 = iS1 + iE1_1 + iE1_2 + iI1 + iA1 + iR1 + iX1;
  double N2 = iS2 + iE2_1 + iE2_2 + iI2 + iA2 + iR2 + iX2;
  
  double patch1 = tau*beta1*( I1*cm[0]+I2*cm[2] )/( N1*cm[0]+N2*cm[2] );
  double patch2 = tau*beta2*( I1*cm[1]+I2*cm[3] )/( N1*cm[1]+N2*cm[3] );
  
  double new_inf1 = iS1*( cm[0]*patch1 + cm[1]*patch2 );
  double new_inf2 = iS2*( cm[2]*patch1 + cm[3]*patch2 );

  // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
  // Rprintf("the value of contact: %f, %f,  %f, %f \n", cm[0], cm[1],cm[2], cm[3] );
   
  double maxtrans1 = R::rpois( new_inf1 );
  double maxtrans2 = R::rpois( new_inf2 );
  double newcase1 = std::min( iS1, maxtrans1 );
  double newcase2 = std::min( iS2, maxtrans2 );
  
  double E1_1toE1_2 = std::min( iE1_1, R::rpois(delta1*tau*iE1_1) ); // 
  double E1_2toI1 = std::min( iE1_2, R::rpois(delta2*tau*iE1_2) ); //
  double I1toR1 = std::min( iI1, R::rpois(gamma*tau*iI1) ); //
  double I1toX1 = std::min( iI1-I1toR1, R::rpois(rate_x*tau*iI1) ); //
  double A1toR1 = std::min( iA1, R::rpois(gamma*tau*iA1) ); //
  double A1toX1 = std::min( iA1-A1toR1, R::rpois(rrxa*rate_x*tau*iA1) ); //
  
  double E2_1toE2_2 = std::min( iE2_1, R::rpois(delta1*tau*iE2_1) ); // 
  double E2_2toI2 = std::min( iE2_2, R::rpois(delta2*tau*iE2_2) ); // 
  double I2toR2 = std::min( iI2, R::rpois(gamma*tau*iI2) ); //
  double A2toR2 = std::min( iA2, R::rpois(gamma*tau*iA2) ); // 
  double I2toX2 = std::min( iI2-I2toR2, R::rpois(rate_x*tau*iI2) ); //
  double A2toX2 = std::min( iA2-A2toR2, R::rpois(rrxa*rate_x*tau*iA2) ); //
  
  // Calculate the change in each state variable
  double dS1 = - newcase1;
  double dE1_1 = + newcase1 - E1_1toE1_2;
  double dE1_2 = E1_1toE1_2 - E1_2toI1;
  double dI1 = (1-fa)*E1_2toI1 - I1toR1 - I1toX1;
  double dA1 = fa*E1_2toI1 - A1toR1 - A1toX1;
  double dR1 = I1toR1 + A1toR1;
  double dX1 = I1toX1 + A1toX1;
  
  double dS2 = - newcase2;
  double dE2_1 = + newcase2 - E2_1toE2_2;
  double dE2_2 = E2_1toE2_2 - E2_2toI2;
  double dI2 = (1-fa)*E2_2toI2 - I2toR2 - I2toX2;
  double dA2 = fa*E2_2toI2 - A2toR2 - A2toX2;
  double dR2 = I2toR2 + A2toR2;
  double dX2 = I2toX2 + A2toX2;
  
  // Update next timestep
  // high-risk population
  SS1[istep+1] = iS1 + dS1;
  EE1_1[istep+1] = iE1_1 + dE1_1;
  EE1_2[istep+1] = iE1_2 + dE1_2;
  II1[istep+1] = iI1 + dI1;
  AA1[istep+1] = iA1 + dA1;
  RR1[istep+1] = iR1 + dR1;
  XX1[istep+1] = iX1 + dX1;
  
  SS2[istep+1] = iS2 + dS2;
  EE2_1[istep+1] = iE2_1 + dE2_1;
  EE2_2[istep+1] = iE2_2 + dE2_2;
  II2[istep+1] = iI2 + dI2;
  AA2[istep+1] = iA2 + dA2;
  RR2[istep+1] = iR2 + dR2;
  XX2[istep+1] = iX2 + dX2;
  
  // cumulative incidence
  CI1[istep+1] = iCI1 + newcase1;
  CI2[istep+1] = iCI2 + newcase2;
  // time in fractional years (ie units parameters are given in)
  time[istep+1] = (istep+1)*tau;
}
// Return results as data.frame
DataFrame sim = DataFrame::create(
  Named("time") = time,
  Named("S1") = SS1,
  Named("E1") = EE1_1 + EE1_2,
  Named("I1") = II1,
  Named("A1") = AA1,
  Named("R1") = RR1,
  Named("X1") = XX1,
  Named("CI1") = CI1,
  
  Named("S2") = SS2,
  Named("E2") = EE2_1 + EE2_2,
  Named("I2") = II2,
  Named("A2") = AA2,
  Named("R2") = RR2,
  Named("X2") = XX2,
  Named("CI2") = CI2
);
return sim;
};

