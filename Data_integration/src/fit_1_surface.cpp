
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






// /** Precision matrix for the anisotropic case, eqn (20) in Lindgren et al. (2011) */    
// namespace R_inla_generalized {
// using namespace Eigen;
// using namespace tmbutils;
// using namespace R_inla;
// 
// template<class Type>
//   SparseMatrix<Type> Q_spde_generalized(spde_t<Type> spde, Type kappa, int alpha=2){
//   Type kappa_pow2 = kappa*kappa;
//   Type kappa_pow4 = kappa_pow2*kappa_pow2;
//   	
//   if( alpha==1 ) return kappa_pow2*spde.M0 + spde.M1;
//   if( alpha==2 ) return kappa_pow4*spde.M0 + Type(2.0)*kappa_pow2*spde.M1 + spde.M2;
// }
// 
// template<class Type>
//   SparseMatrix<Type> Q_spde_generalized(spde_aniso_t<Type> spde, Type kappa, matrix<Type> H, int alpha=2){
// 
//   int i;
//   Type kappa_pow2 = kappa*kappa;
//   Type kappa_pow4 = kappa_pow2*kappa_pow2;
//   
//   int n_s = spde.n_s;
//   int n_tri = spde.n_tri;
//   vector<Type> Tri_Area = spde.Tri_Area;
//   matrix<Type> E0 = spde.E0;
//   matrix<Type> E1 = spde.E1;
//   matrix<Type> E2 = spde.E2;
//   matrix<int> TV = spde.TV;
//   SparseMatrix<Type> G0 = spde.G0;
//   SparseMatrix<Type> G0_inv = spde.G0_inv;
// 	  	  
//   //Type H_trace = H(0,0)+H(1,1);
//   //Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
//   SparseMatrix<Type> G1_aniso(n_s,n_s); 
//   SparseMatrix<Type> G2_aniso(n_s,n_s); 
//   // Calculate adjugate of H
//   matrix<Type> adj_H(2,2);
//   adj_H(0,0) = H(1,1);
//   adj_H(0,1) = -1 * H(0,1);
//   adj_H(1,0) = -1 * H(1,0);
//   adj_H(1,1) = H(0,0);
//   // Calculate new SPDE matrices
// 
//   // Calculate G1 - pt. 1
//   array<Type> Gtmp(n_tri,3,3);
//   for(i=0; i<n_tri; i++){    
//     // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
//     Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
//     Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
//     Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
//     Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
//     Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
//     Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
//   }
//   // Calculate G1 - pt. 2
//   for(i=0; i<n_tri; i++){
//     G1_aniso.coeffRef(TV(i,1),TV(i,0)) = G1_aniso.coeffRef(TV(i,1),TV(i,0)) + (Gtmp(i,0,1));  
//     G1_aniso.coeffRef(TV(i,0),TV(i,1)) = G1_aniso.coeffRef(TV(i,0),TV(i,1)) + (Gtmp(i,0,1));  
//     G1_aniso.coeffRef(TV(i,2),TV(i,1)) = G1_aniso.coeffRef(TV(i,2),TV(i,1)) + (Gtmp(i,1,2));  
//     G1_aniso.coeffRef(TV(i,1),TV(i,2)) = G1_aniso.coeffRef(TV(i,1),TV(i,2)) + (Gtmp(i,1,2));  
//     G1_aniso.coeffRef(TV(i,2),TV(i,0)) = G1_aniso.coeffRef(TV(i,2),TV(i,0)) + (Gtmp(i,0,2));  
//     G1_aniso.coeffRef(TV(i,0),TV(i,2)) = G1_aniso.coeffRef(TV(i,0),TV(i,2)) + (Gtmp(i,0,2));  
//     G1_aniso.coeffRef(TV(i,0),TV(i,0)) = G1_aniso.coeffRef(TV(i,0),TV(i,0)) + (Gtmp(i,0,0));  
//     G1_aniso.coeffRef(TV(i,1),TV(i,1)) = G1_aniso.coeffRef(TV(i,1),TV(i,1)) + (Gtmp(i,1,1));  
//     G1_aniso.coeffRef(TV(i,2),TV(i,2)) = G1_aniso.coeffRef(TV(i,2),TV(i,2)) + (Gtmp(i,2,2));  
//   }
//   G2_aniso = G1_aniso * G0_inv * G1_aniso; 
// 
//   if( alpha==1 ) return kappa_pow2*G0 + G1_aniso;
//   if( alpha==2 ) return kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso;
// }
// } // end namespace R_inla

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
  DATA_VECTOR( C_i );       	// Vector of responses (counts) of each species at each sampled location 
  DATA_VECTOR( P_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( A_s );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
  DATA_IVECTOR( S_i ); // Site index for each sample
  DATA_MATRIX( X_s );  //design matrix for fixed effects
  
  // normalization flag - used for speed -up
  DATA_INTEGER(flag); // flag == 0 => no data contribution added to jnll
  
  // Indices
  DATA_INTEGER( n_i ); // Number of data points in space
  DATA_INTEGER( n_s ); // Number of grid cells
  DATA_INTEGER( n_eta ); //mesh points in space mesh
  

  // SPDE objects
  DATA_SPARSE_MATRIX ( M0 );
  DATA_SPARSE_MATRIX ( M1 );
  DATA_SPARSE_MATRIX ( M2 );
  //DATA_SPARSE_MATRIX ( Aproj );
  DATA_IVECTOR ( Eta_index );
  
  // Options
  DATA_VECTOR ( options );
  // options [0] == 1 : use normalization trick

  // Prior specifications
  DATA_VECTOR( beta_pri );    //NOT CURRENTLY USED
  DATA_VECTOR( matern_pri );  //NOT CURRENTLY USED
  // matern_pri = c(a, b, c, d): P( range < a) = b; P( sigma > c) = d
  Type matern_par_a = matern_pri[0]; // range limit : rho0
  Type matern_par_b = matern_pri[1]; // range prob : alpha_rho
  Type matern_par_c = matern_pri[2]; // field sd limit : sigma0
  Type matern_par_d = matern_pri[3]; // field sd prob : alpha_sigma

 
  // Parameters 
  PARAMETER_VECTOR( Beta );              // fixed effects on density
  PARAMETER( log_tau );      //spatial precision
  PARAMETER( log_kappa );
  PARAMETER_VECTOR( Etainput_s );           // spatial random effects - (n_sp x n_eta)

    // derived sizes
  int n_b = X_s.row(0).size();
  
  // global stuff
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  Type jnll = 0.0;
  
  // Make spatial precision matrix
  SparseMatrix <Type > Q_ss = spde_Q( log_kappa , log_tau , M0 , M1 , M2);
  
  // Transform some of our parameters
  Type sp_range = sqrt(8.0) / exp(log_kappa);
  Type sp_sigma = 1.0 / sqrt (4.0 * 3.14159265359 * exp(2.0*log_tau) * exp(2.0*log_kappa));
  

  // Transform random effects  #not sure we need this anymore; Osgood-Zimmerman & Wakefield don't have it in there
  // vector<Type> Eta_s(n_eta);
  // for(int ieta=0;ieta<n_eta;ieta++){
  //   Eta_s(ieta) = Etainput_s(ieta) / exp(logtau_eta);
  // }

  // the random effects . we do this first so to do the
  // normalization outside of every optimization step
  // NOTE : likelihoods from namespace ’density ’ already return NEGATIVE
  // log - liks so we add other likelihoods return positive log - liks
  if( options [0] == 1){
    // then we are not calculating the normalizing constant in the inner opt
    // that norm constant means taking an expensive determinant of Q_ss
    jnll_comp(0) += GMRF(Q_ss , false )(Etainput_s);
    // return without data ll contrib to avoid unneccesary log ( det (Q)) calcs
    if ( flag == 0) return jnll_comp(0);
  } 
  else {
    jnll_comp(0) += GMRF(Q_ss)(Etainput_s);
  }
  
  
  // Predicted densities
  vector<Type> Z_s(n_s);
  vector<Type> log_Z_s(n_s);
  vector<Type> E_count_obs(n_i);
  vector<Type> RE(n_s);
  vector<Type> linpredZ_s(n_s);
  linpredZ_s = X_s * Beta;

  for(int is=0;is<n_s;is++){
    RE(is)=Etainput_s(Eta_index(is));
    log_Z_s(is) = linpredZ_s(is) + RE(is);
    Z_s(is) = exp( log_Z_s(is));
  }
  

  // Probability of counts
  E_count_obs = E_count_obs.setZero();
  for(int i=0; i<n_i; i++){
    E_count_obs(i)=P_i(i)*Z_s(S_i(i));
  }
  
  for(int i=0;i<n_i;i++){
      jnll_comp(1) -= dpois( C_i(i), E_count_obs(i), true );
  }
  
  //priors
  //jnll_comp(2) -= dnorm (Beta , beta_pri [0] , beta_pri [1] , true );
  

  // Total objective
   jnll = jnll_comp.sum();

  //std::cout<<"Range "<<Range_eta<<"\n";
   REPORT( Z_s );
   REPORT( Beta );
   REPORT( jnll_comp );
   REPORT( jnll );
   REPORT( E_count_obs);
   REPORT( sp_sigma );
   REPORT( sp_range );
   REPORT( Etainput_s );
   REPORT( log_tau );
   REPORT( log_kappa );
   REPORT( RE );
   REPORT( Z_s );

  // Bias correction output
  ADREPORT(log_Z_s);

  return jnll;
}
