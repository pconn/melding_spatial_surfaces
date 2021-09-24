# integrate_SSL.R  integrate POP and UDs for SSL

library(TMB)
library(Matrix)

load("SSL_TMB_data.RData")
TmbFile = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_multiple_surfaces_work"
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 
dyn.load( dynlib(TmbFile) )

Data_SSL$Omega_Y1 = as(Data_SSL$Omega_Y[,,1]*diag(Data_SSL$n_s),"dgTMatrix")
Data_SSL$Omega_Y2 = as(Data_SSL$Omega_Y[,,2]*diag(Data_SSL$n_s),"dgTMatrix")
which_el = which(names(Data_SSL)=="Omega_Y")
Data_SSL[[which_el]]<-NULL
gc()  #clean up workspace

Params_all <- list("log_tau_mu"=0, "log_kappa_mu"=0,"Mu_s"=rep(0,Data_SSL$n_mesh))      
Params_all$Xi_s = matrix(0,Data_SSL$n_mesh,2)
Params_all$log_alpha = Params_all$log_tau_xi = Params_all$log_kappa_xi = rep(0,2)
#Random_all <- c("Mu_s","Xi_s")
Random_all <- "Mu_s"
Map_all = list("Xi_s"=factor(matrix(NA,Data_SSL$n_mesh,2)),"log_kappa_xi"=factor(rep(NA,2)),
               "log_tau_xi"=factor(rep(NA,2)))
#Map_all=NULL
#n_surf=2
#Map_all = list("Xi_s"=factor(c(rep(NA,Data_SSL$n_mesh),c((2*n_surf-1):(Data_SSL$n_mesh*(n_surf-1)+2*n_surf-2)))),"log_kappa_xi"=factor(c(NA,c(1:(n_surf-1)))),
#               "log_tau_xi"=factor(c(NA,c(n_surf:(2*n_surf-2)))))

Obj = MakeADFun( data=Data_SSL, parameters=Params_all, random=Random_all, map=Map_all, DLL="fit_multiple_surfaces_work",silent=FALSE)
st_time = Sys.time()
Obj$fn( Obj$par )
cat(Sys.time()-st_time)
#Obj <- normalize ( Obj , flag ="flag", value = 0)
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
Report = Obj$report()
