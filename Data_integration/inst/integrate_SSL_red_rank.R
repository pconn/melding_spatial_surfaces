# integrate_SSL.R  integrate POP and UDs for SSL

library(TMB)
library(Matrix)

load("SSL_TMB_data.RData")
TmbFile = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_multiple_surfaces_rr"
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 
dyn.load( dynlib(TmbFile) )

Params_all <- list("log_tau_mu"=0, "log_kappa_mu"=0,"Mu_s"=rep(0,Data_SSL$n_mesh))      
Params_all$Xi_s = matrix(0,Data_SSL$n_mesh,2)
Params_all$log_alpha = Params_all$log_tau_xi = Params_all$log_kappa_xi = rep(0,2)

### 1) RUN model with process errors (xi) set to zero
Random_all <- "Mu_s"
Map_all = list("Xi_s"=factor(matrix(NA,Data_SSL$n_mesh,2)),"log_kappa_xi"=factor(rep(NA,2)),
               "log_tau_xi"=factor(rep(NA,2)))

Obj = MakeADFun( data=Data_SSL, parameters=Params_all, random=Random_all, map=Map_all, DLL="fit_multiple_surfaces_rr",silent=FALSE)
Obj$fn( Obj$par )
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
Report = Obj$report()

save.image("output_rr_nobias.RData")

### 2) RUN model with process errors (xi) estimated

Map_all=NULL
Random_all <- c("Mu_s","Xi_s")

Obj = MakeADFun( data=Data_SSL, parameters=Params_all, random=Random_all, map=Map_all, DLL="fit_multiple_surfaces_rr",silent=FALSE)
Obj$fn( Obj$par )
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
Report = Obj$report()

save.image("output_rr_bias.RData")


