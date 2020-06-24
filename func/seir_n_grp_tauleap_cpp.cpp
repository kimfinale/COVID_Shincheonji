#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


List seir_3grp_tauleap(List params) {

// chained operations are tricky in cpp
// pull out list w/in list into its own object
List init = params["init"];

// use Rcpp as() function to "cast" R vector to cpp scalar
int nsteps = as<int>(params["nsteps"]);

// ints (IntegerVector), since rpois returns double,

NumericVector SS1( nsteps, init["S1"] );
NumericVector EE1_1( nsteps, init["E1_1"] );
NumericVector EE1_2( nsteps, init["E1_2"] );
NumericVector II1( nsteps, init["I1"] );
NumericVector AA1( nsteps, init["A1"] );//asymptomatic
NumericVector XX1( nsteps, init["X1"] );//isolated
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


double ngroups = params["ngroups"];
NumericVector betas = params[ "betas" ]; // cumulative

NumericMatrix pops_S( nsteps, ngroups ); 
NumericMatrix pops_E1( nsteps, ngroups ); 
NumericMatrix pops_E2( nsteps, ngroups ); 
NumericMatrix pops_I( nsteps, ngroups ); 
NumericMatrix pops_R( nsteps, ngroups ); 
NumericMatrix pops_X( nsteps, ngroups ); 
// fill time w/zeros
NumericVector time(nsteps);

// pull out params for easy reading 
NumericVector cm = params[ "cm" ];//contact matrix
double delta = params["delta"];
double gamma = params["gamma"];
double tau = params["tau"];
double fa = params["frac_asymp"];
// double rate_x = params["rate_x"];
double rrx1a = params["rel_rate_isol1_asymp"]; //relative isolation rate for asymptomatic
double rrx2a = params["rel_rate_isol2_asymp"]; //relative isolation rate for asymptomatic
double rrx3a = rrx2a;
double time_intervention_start = params["time_intervention_start"];
double daily_frac_reduc1 = params[ "daily_frac_reduc1" ];
double daily_frac_reduc2 = params[ "daily_frac_reduc2" ];

double min_delay1 = params[ "min_delay1" ];
double min_delay2 = params[ "min_delay2" ];
double min_beta1 = params[ "min_beta1" ];
double min_beta2 = params[ "min_beta2" ];

double pop3_contact_start = params[ "pop3_contact_start" ];
double base_behavior_start = params[ "base_behavior_start" ];
// Rprintf("the value of beta[1] : %f \n", beta1);
// Rprintf("the value of beta[2] : %f \n", beta2);

// Calculate the number of events for each step, update state vectors
for( int istep = 0; istep < (nsteps-1); istep++ ) {
  // pull out this step's scalars for easier reading
  // and to avoid compiler headaches
  double iS1 = SS1[istep]; //  rpois returns double
  double iE1_1 = EE1_1[istep];
  double iE1_2 = EE1_2[istep];
  double iI1 = II1[istep];
  double iA1 = AA1[istep];
  double iR1 = RR1[istep];
  double iX1 = XX1[istep];
  double iCI1 = CI1[istep];
  
  double iS2 = SS2[istep];
  double iE2_1 = EE2_1[istep];
  double iE2_2 = EE2_2[istep];
  double iI2 = II2[istep];
  double iA2 = AA2[istep];
  double iX2 = XX2[istep];
  double iR2 = RR2[istep];
  double iCI2 = CI2[istep];
  
  NumericVector step_S = pop_S( istep, _ );
  NumericVector step_E1 = pop_E1( istep, _ );
  NumericVector step_E2 = pop_E2( istep, _ );
  NumericVector step_I = pop_I( istep, _ );
  NumericVector step_R = pop_R( istep, _ );
  NumericVector step_X = pop_X( istep, _ );
  
  // for( int g = 0; g < ngroups; g++ ){
  //   npops[ istep, g ] = 
  // }
  double iS3 = SS3[ istep ];
  double iE3_1 = EE3_1[ istep ];
  double iE3_2 = EE3_2[ istep ];
  double iI3 = II3[ istep ];
  double iA3 = AA3[ istep ];
  double iX3 = XX3[ istep ];
  double iR3 = RR3[ istep ];
  double iCI3 = CI3[ istep ];
  
  double beta1 = params[ "beta1" ];
  double beta2 = params[ "beta2" ];
  double beta3 = params[ "beta3" ];
  
  double rate_x1 = params["rate_isol1"]; //initialize
  double delay1 = 1/rate_x1;
  double rate_x2 = params[ "rate_isol2" ];
  double delay2 = 1/rate_x2;
  double rate_x3 = rate_x2;
  /////////////////////////
  // State Equations
  /////////////////////////
  
  // R::rpois always returns a single value
  // to return multiple (e.g. Integer/NumericVector, 
  // use Rcpp::rpois(int ndraw, param) and friends
  // Prevent negative states

  double today = (istep+1)*tau;
  double time_diff = today - time_intervention_start;
  
  daily_frac_reduc2 = daily_frac_reduc1; // special enforcing in high risk group?
  
  // if( time_diff > 0 && today < pop3_contact_start ){
  if( time_diff > 0 && today < base_behavior_start ){
    delay1 = std::max( delay1*(1 - daily_frac_reduc1*time_diff), min_delay1 );
    delay2 = std::max( delay2*(1 - daily_frac_reduc2*time_diff), min_delay2 );
    beta1 = std::max( beta1*(1 - daily_frac_reduc1*time_diff), min_beta1 );
    beta2 = std::max( beta2*(1 - daily_frac_reduc2*time_diff), min_beta2 );
  }
  
  if( today < pop3_contact_start ){
    cm[2] = 0;
    cm[5] = 0;
  	cm[6] = 0;
  	cm[7] = 0;
  	cm[8] = 0;
  }
  
  rate_x1 = 1/delay1;
  rate_x2 = 1/delay2;
  rate_x3 = 1/delay2;
  
  double I1 = iI1 + iA1; // asymptomatic individuals have the same transmissibility
  double I2 = iI2 + iA2;
  double I3 = iI3 + iA3;
  
  double N1 = iS1 + iE1_1 + iE1_2 + iI1 + iA1 + iR1 + iX1;
  double N2 = iS2 + iE2_1 + iE2_2 + iI2 + iA2 + iR2 + iX2;
  double N3 = iS3 + iE3_1 + iE3_2 + iI3 + iA3 + iR3 + iX3;
  
  double cr1 = (I1*cm[0]+I2*cm[3]+I3*cm[6])/(N1*cm[0]+N2*cm[3]+N3*cm[6]);
  double cr2 = (I1*cm[1]+I2*cm[4]+I3*cm[7])/(N1*cm[1]+N2*cm[4]+N3*cm[7]); 
  double cr3 = 0;
  if( (cm[2] + cm[5] + cm[8]) > 0  ){
    cr3 = (I1*cm[2]+I2*cm[5]+I3*cm[8])/(N1*cm[2]+N2*cm[5]+N3*cm[8]);
  }
  double foi1 = tau*(beta1*cm[0]*cr1 + beta2*cm[1]*cr2 + beta3*cm[2]*cr3);
  double foi2 = tau*(beta1*cm[3]*cr1 + beta2*cm[4]*cr2 + beta3*cm[5]*cr3);
  double foi3 = tau*(beta1*cm[6]*cr1 + beta2*cm[7]*cr2 + beta3*cm[8]*cr3);
  
  double new_inf1 = iS1*foi1;
  double new_inf2 = iS2*foi2;
  double new_inf3 = iS3*foi3;
  
  // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
  // Rprintf("the value of contact: %f, %f,  %f, %f \n", cm[0], cm[1],cm[2], cm[3] );
   
  double maxtrans1 = R::rpois( new_inf1 );
  double maxtrans2 = R::rpois( new_inf2 );
  double maxtrans3 = R::rpois( new_inf3 );
  
  double newcase1 = std::min( iS1, maxtrans1 );
  double newcase2 = std::min( iS2, maxtrans2 );
  double newcase3 = std::min( iS3, maxtrans3 );
  
  double E1_1toE1_2 = std::min( iE1_1, R::rpois(2*delta*tau*iE1_1) ); // 
  double E1_2toI1 = std::min( iE1_2, R::rpois(2*delta*tau*iE1_2) ); // 
  double I1toR1 = std::min( iI1, R::rpois(gamma*tau*iI1) ); //
  double I1toX1 = std::min( iI1-I1toR1, R::rpois(rate_x1*tau*iI1) ); //
  double A1toR1 = std::min( iA1, R::rpois(gamma*tau*iA1) ); //
  double A1toX1 = std::min( iA1-A1toR1, R::rpois(rrx1a*rate_x1*tau*iA1) ); //
  
  double E2_1toE2_2 = std::min( iE2_1, R::rpois(2*delta*tau*iE2_1) ); // 
  double E2_2toI2 = std::min( iE2_2, R::rpois(2*delta*tau*iE2_2) ); // 
  double I2toR2 = std::min( iI2, R::rpois(gamma*tau*iI2) ); //
  double A2toR2 = std::min( iA2, R::rpois(gamma*tau*iA2) ); // 
  double I2toX2 = std::min( iI2-I2toR2, R::rpois(rate_x2*tau*iI2) ); //
  double A2toX2 = std::min( iA2-A2toR2, R::rpois(rrx2a*rate_x2*tau*iA2) ); //

  double E3_1toE3_2 = std::min( iE3_1, R::rpois(2*delta*tau*iE3_1) ); // 
  double E3_2toI3 = std::min( iE3_2, R::rpois(2*delta*tau*iE3_2) ); // 
  double I3toR3 = std::min( iI2, R::rpois(gamma*tau*iI3) ); //
  double A3toR3 = std::min( iA2, R::rpois(gamma*tau*iA3) ); // 
  double I3toX3 = std::min( iI2-I2toR2, R::rpois(rate_x3*tau*iI3) ); //
  double A3toX3 = std::min( iA2-A2toR2, R::rpois(rrx3a*rate_x3*tau*iA2) ); // 
  
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
 
  double dS3 = - newcase3;
  double dE3_1 = + newcase3 - E3_1toE3_2;
  double dE3_2 = E3_1toE3_2 - E3_2toI3;
  double dI3 = (1-fa)*E3_2toI3 - I3toR3 - I3toX3;
  double dA3 = fa*E3_2toI3 - A3toR3 - A3toX3;
  double dR3 = I3toR3 + A3toR3;
  double dX3 = I3toX3 + A3toX3;
  
  // Update next timestep
  // high-risk population
  SS1[istep+1] = iS1 + dS1;
  EE1_1[istep+1] = iE1_1 + dE1_1;
  EE1_2[istep+1] = iE1_2 + dE1_2;
  II1[istep+1] = iI1 + dI1;
  AA1[istep+1] = iA1 + dA1;
  RR1[istep+1] = iR1 + dR1;
  XX1[istep+1] = iX1 + dX1;
  CI1[istep+1] = iCI1 + newcase1;// cumulative incidence
  
  SS2[istep+1] = iS2 + dS2;
  EE2_1[istep+1] = iE2_1 + dE2_1;
  EE2_2[istep+1] = iE2_2 + dE2_2;
  II2[istep+1] = iI2 + dI2;
  AA2[istep+1] = iA2 + dA2;
  RR2[istep+1] = iR2 + dR2;
  XX2[istep+1] = iX2 + dX2;
  CI2[istep+1] = iCI2 + newcase2;
  
  SS3[istep+1] = iS3 + dS3;
  EE3_1[istep+1] = iE3_1 + dE3_1;
  EE3_2[istep+1] = iE3_2 + dE3_2;
  II3[istep+1] = iI3 + dI3;
  AA3[istep+1] = iA3 + dA3;
  RR3[istep+1] = iR3 + dR3;
  XX3[istep+1] = iX3 + dX3;
  CI3[istep+1] = iCI3 + newcase3;
  
  // time in fractional years (ie units parameters are given in)
  time[istep+1] = (istep+1)*tau;
}

DataFrame sim = DataFrame::create(
  Named("time") = time,
  Named("S1") = SS1,
  Named("E1") = EE1_1 + EE1_2,
  Named("I1") = II1,
  Named("A1") = AA1,
  Named("R1") = RR1,
  Named("X1") = XX1,
  Named("CI1") = CI1,
  
  Named("S2") = SS2 + SS3,
  Named("E2") = EE2_1 + EE2_2 + EE3_1 + EE3_2,
  Named("I2") = II2 + II3,
  Named("A2") = AA2 + AA3,
  Named("R2") = RR2 + RR3,
  Named("X2") = XX2 + XX3,
  Named("CI2") = CI2 + CI3
  // 
  // Named("S3") = SS3,
  // Named("E3") = EE3_1 + EE3_2,
  // Named("I3") = II3,
  // Named("A3") = AA3,
  // Named("R3") = RR3,
  // Named("X3") = XX3,
  // Named("CI3") = CI3
);
return sim;
};


