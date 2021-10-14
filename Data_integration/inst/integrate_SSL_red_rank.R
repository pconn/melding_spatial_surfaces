# integrate_SSL.R  integrate POP and UDs for SSL

library(TMB)
library(Matrix)

load("SSL_TMB_data.RData")
TmbFile = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_multiple_surfaces_rr"
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 
dyn.load( dynlib(TmbFile) )

Data_SSL$Omega_Y1 = as(Data_SSL$Omega_Y[,,1]*diag(Data_SSL$n_s),"dgTMatrix")
Data_SSL$Omega_Y2 = as(Data_SSL$Omega_Y[,,2]*diag(Data_SSL$n_s),"dgTMatrix")
#Data_SSL$Omega_Y1 = as(Data_SSL$Omega_Y[,,1],"dgTMatrix")
#Data_SSL$Omega_Y2 = as(Data_SSL$Omega_Y[,,2],"dgTMatrix")


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
n_surf=2
#Map_all = list("Xi_s"=factor(c(rep(NA,Data_SSL$n_mesh),c((2*n_surf-1):(Data_SSL$n_mesh*(n_surf-1)+2*n_surf-2)))),"log_kappa_xi"=factor(c(NA,c(1:(n_surf-1)))),
#               "log_tau_xi"=factor(c(NA,c(n_surf:(2*n_surf-2)))))

#Data_SSL2 = list(flag=Data_SSL$flag,n_s=Data_SSL$n_s,n_mesh=Data_SSL$n_mesh,n_surf=Data_SSL$n_surf,options=Data_SSL$options,matern_pri=Data_SSL$matern_pri,
#                 M0=Data_SSL$M0,M1=Data_SSL$M1,M2=Data_SSL$M2,Y_i=Data_SSL$Y_i,Omega_Y1=Data_SSL$Omega_Y1,Omega_Y2=Data_SSL$Omega_Y2,A=Data_SSL$A)

Obj = MakeADFun( data=Data_SSL, parameters=Params_all, random=Random_all, map=Map_all, DLL="fit_multiple_surfaces_rr",silent=FALSE)
st_time = Sys.time()
Obj$fn( Obj$par )
cat(Sys.time()-st_time)
#Obj <- normalize ( Obj , flag ="flag", value = 0)
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
Report = Obj$report()

save.image("output_rr_nobias.RData")

load("output_rr_bias.RData")
#load("output_rr_nobias.RData")
#plot surface
load("./data_integration/data/POPfits.RDa")
Grid_locs = POPfits[,c(1:2)]
source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')
plot_N_map_xy(N=Data_SSL$Y_i[,1],XY=Grid_locs,leg.title="POP")
plot_N_map_xy(N=Data_SSL$Y_i[,2],XY=Grid_locs,leg.title="UD")
plot_N_map_xy(N=Report$Mu_pred,XY=Grid_locs,leg.title="Mu_s")

POP_abs = exp(Data_SSL$Y_i[,1])/sum(exp(Data_SSL$Y_i[,1]))
UD_abs = exp(Data_SSL$Y_i[,2])/sum(exp(Data_SSL$Y_i[,2]))
Mean_abs = 0.5*(POP_abs+UD_abs)
plot_N_map_xy(N=POP_abs,XY=Grid_locs,leg.title="POP-abs")
plot_N_map_xy(N=UD_abs,XY=Grid_locs,leg.title="UD-abs")
plot_N_map_xy(N=Report$Pi,XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=Mean_abs,XY=Grid_locs,leg.title="Mean-abs")

#plot_N_map_xy(N=Report$Mu_s[1:13658],XY=Grid_locs,leg.title="Mu_s")
plot_N_map_xy(N=Report$Xi_pred[,1],XY=Grid_locs,leg.title="Bias-POP")
plot_N_map_xy(N=Report$Xi_pred[,2],XY=Grid_locs,leg.title="Bias-UD")


plot_N_map_xy(N=exp(Data_SSL$Y_i[,1])/sum(exp(Data_SSL$Y_i[,1])),XY=Grid_locs,leg.title="POP")
plot_N_map_xy(N=exp(Data_SSL$Y_i[,2])/sum(exp(Data_SSL$Y_i[,2])),XY=Grid_locs,leg.title="UD")
plot_N_map_xy(N=Report$Pi,XY=Grid_locs,leg.title="Rel abundance")


