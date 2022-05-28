
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function to make sparse SPDE precision matrix
// Inputs :
// logkappa : log ( kappa ) parameter value
// logtau : log ( tau ) parameter value
// M0 , M1 , M2: these sparse matrices are output from :
// R:: INLA :: inla . spde2 . matern () $ param . inla $M*
template <class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0 ,
                              SparseMatrix <Type> M1 , SparseMatrix <Type> M2) {
  SparseMatrix <Type> Q;
  Type kappa2 = exp (2. * logkappa);
  Type kappa4 = kappa2 * kappa2;
  Q = pow (exp(logtau), 2.) * (kappa4*M0 + Type(2.0)*kappa2*M1+M2);
  return Q;
}

// helper function to use the same penalized complexity prior on
// matern params that is used in INLA

template <class Type>
Type dPCPriSPDE (Type logtau, Type logkappa,
                  Type matern_par_a, Type matern_par_b,
                  Type matern_par_c, Type matern_par_d,
                  int give_log=0)
{

  Type penalty; // prior contribution to jnll
  
  Type d=2.; // dimension
  Type lambda1= -log(matern_par_b)*pow(matern_par_a, d/2.);
  Type lambda2= -log(matern_par_d)/matern_par_c;
  Type range = sqrt(8.0)/exp(logkappa);
  Type sigma = 1.0/sqrt(4.0*3.14159265359*exp(2.0*logtau)*
    exp(2.0*logkappa));
  
  penalty=(-d/2. - 1.)*log(range)-lambda1*pow(range,-d/2.)-lambda2*sigma;
  // Note : (rho , sigma ) --> (x= log kappa , y= log tau ) -->
  // transforms : rho = sqrt (8) /e^x & sigma = 1/( sqrt (4 pi)*e^x*e^y)
  // --> Jacobian : |J| propto e^( -y -2x)
  Type jacobian = -logtau-2.0*logkappa;
  penalty += jacobian;
  
  if( give_log ) return penalty; else return exp(penalty);
}



//dmvnorm function with pre-calculated matrix inverse, log determinant - ignore normalizing constant
template <class Type>
Type neg_log_dmvnorm(vector<Type> x, vector<Type> mu, matrix<Type> omega) {
 vector<Type> ss = x-mu;
 return ((0.5*ss*(omega*ss)).sum());
}

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}


template<class Type>
Type objective_function<Type>::operator() ()
{


  // Data
  DATA_MATRIX( Y_i );       	// (max(n_s) * n_surf) MATRIX of surface values for each data type at each sampled location 
                              // each column holds surface values for surface i; these can be ragged
  DATA_ARRAY(Omega_Y);    //inverse of Sigma_Y (log scale).  Dimensions are assumed to be max(n_s),max(n_s),n_surf
  DATA_ARRAY(A_i);        // array holding weighting matrices for converting mu into y; dimensions are (max(n_s),n_s_mu,n_surf)

  // normalization flag - used for speed -up
  DATA_INTEGER(flag); // flag == 0 => no data contribution added to jnll
  
  // Indices
  DATA_INTEGER( n_s_mu); //number of grid cells in Mu
  DATA_IVECTOR( n_s ); // Number of grid cells in each relative abundance surface
  DATA_INTEGER( n_surf );
  DATA_IVECTOR( n_knots_xi ); //# points in space mesh for xi (set to zero if xi not estimated for a given surface)
  
  // SPDE objects  - assume these are all the same for each surface
  DATA_SPARSE_MATRIX ( M0_mu );
  DATA_SPARSE_MATRIX ( M1_mu );
  DATA_SPARSE_MATRIX ( M2_mu );
  DATA_SPARSE_MATRIX ( M0 );
  DATA_SPARSE_MATRIX ( M1 );
  DATA_SPARSE_MATRIX ( M2 );
  
  // Options
  DATA_VECTOR ( options );
  // options [0] == 1 : use normalization trick



 
  // Parameters 
  PARAMETER_VECTOR( log_alpha); //scaling params for each surface
  PARAMETER( log_tau_mu );      //spatial precision
  PARAMETER( log_kappa_mu );
  PARAMETER_VECTOR( log_tau_xi );      //spatial precision; one for each surface
  PARAMETER_VECTOR( log_kappa_xi );

  PARAMETER_VECTOR( Mu_s );           // spatial mean RE 
  PARAMETER_MATRIX( Xi_s );           // spatial bias REs 

  // global stuff
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  Type jnll = 0.0;
  vector<Type>alpha=exp(log_alpha);
  int n_s_max = A_i.dim(0);
  
  
  // Make spatial precision matrix
  SparseMatrix <Type > Q_mu = spde_Q( log_kappa_mu , log_tau_mu , M0_mu , M1_mu , M2_mu);
  SparseMatrix <Type > Q_xi; 
  // the random effects . we do this first so to do the
  // normalization outside of every optimization step
  if( options [0] == 1){
    // then we are not calculating the normalizing constant in the inner opt
    // that norm constant means taking an expensive determinant of Q_ss
    jnll_comp(0) += GMRF(Q_mu , false )(Mu_s);
    for(int isurf=0;isurf<n_surf;isurf++){
      if(n_knots_xi(isurf)>0){
        Q_xi = spde_Q( log_kappa_xi(isurf) , log_tau_xi(isurf) , M0 , M1 , M2);
        jnll_comp(0) += GMRF(Q_xi , false )(Xi_s.col(isurf));
      }
    }
    // return without data ll contrib to avoid unneccesary log ( det (Q)) calcs
    if ( flag == 0) return jnll_comp(0) ; //see documentation for TMB::normalize
  } 
  else {
    jnll_comp(0) += GMRF(Q_mu)(Mu_s);
    for(int isurf=0;isurf<n_surf;isurf++){
      if(n_knots_xi(isurf)>0){
        Q_xi = spde_Q( log_kappa_xi(isurf) , log_tau_xi(isurf) , M0 , M1 , M2);
        jnll_comp(0) += GMRF(Q_xi)(Xi_s.col(isurf));
      }
    }
    
  }
  
 
  matrix<Type> Expected(n_s_max,n_surf);
  vector<Type> Expected_tmp(n_s(2));
  vector<Type> Y_i_tmp(n_s(2));
  matrix<Type> Omega_tmp(n_s(2),n_s(2));
  //matrix<Type> Y_i_cur(n_s_max,n_surf);
  //matrix<Type> Expected_cur(n_s_max,n_surf);
  matrix<Type> Omega_cur(n_s_max,n_s_max);
  Type cursum;
  // MVN likelihoods
  for(int isurf=0;isurf<n_surf;isurf++){
    for(int is=0;is<n_s(isurf);is++){
      cursum=0;
      for(int imu=0; imu<n_s_mu;imu++){
        cursum+=A_i(is,imu,isurf)*exp(Mu_s(imu));
      }
      Expected(is,isurf) = alpha(isurf) + log(cursum) + Xi_s(is,isurf);
      
      //Expected_cur(is)=Expected(is,isurf);
      //Y_i_cur(is)=Y_i(is,isurf);
      for(int is2=0;is2<n_s(isurf);is2++){
        Omega_cur(is,is2)=Omega_Y(is,is2,isurf);
      }
      
    }
    vector<Type> Y_i_cur = Y_i.col(isurf).segment(0,n_s(isurf));
    vector<Type> Expected_cur = Expected.col(isurf).segment(0,n_s(isurf));
    matrix<Type> Omega2 = Omega_cur.block(0,0,n_s(isurf),n_s(isurf));
    jnll_comp(1) += neg_log_dmvnorm(Y_i_cur, Expected_cur, Omega2);
    
    Omega_tmp = Omega2;
    Expected_tmp = Expected_cur;
    Y_i_tmp = Y_i_cur;
  }
  

  //priors
  //jnll_comp(2) -= dnorm (Beta , beta_pri [0] , beta_pri [1] , true );

  Type mu_sum = 0;
  for(int is=0;is<n_s_mu;is++)mu_sum+=exp(Mu_s(is));
  vector<Type> Pi(n_s_mu); 
  for(int is=0;is<n_s_mu;is++)Pi(is)= exp(Mu_s(is))/mu_sum;
  
  // Total objective
   jnll = jnll_comp.sum();


   REPORT( jnll_comp );
   REPORT( jnll );
   REPORT( Expected);
   REPORT( Xi_s );
   REPORT( Mu_s );
   REPORT( log_tau_mu );
   REPORT( log_kappa_mu );
   REPORT( log_tau_xi );
   REPORT( log_kappa_xi );
   REPORT( Pi );
   REPORT( alpha );
   REPORT( Omega_tmp);
   REPORT( Expected_tmp);
   REPORT( Y_i_tmp);
 
  // Bias correction output
  ADREPORT(Pi);

  return jnll;
}
