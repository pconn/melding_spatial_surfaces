# simulations for combining spatial surfaces

library( TMB )
library( INLA )
library( mvtnorm)

get_deriv_mat_pi <- function(Y){  
  Y_sum = sum(Y)
  Y_sum2 = Y_sum^2
  D_mat = matrix(0,length(Y),length(Y))
  for(irow in 1:length(Y))D_mat[irow,]=-Y
  diag(D_mat)= Y_sum - Y
  D_mat = D_mat/Y_sum2
  diag(D_mat)=diag(D_mat)
  D_mat
}

TmbFile1 = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_1_surface"
compile(paste0(TmbFile1,".cpp"),"-O1 -g",DLLFLAGS="") 

TmbFile2 = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_multiple_surfaces"
compile(paste0(TmbFile2,".cpp"),"-O1 -g",DLLFLAGS="") 

dyn.load( dynlib(TmbFile1) )
dyn.load( dynlib(TmbFile2) )


init.seed = 2020
N=10000 #abundance

n_cols = 15
n_rows = 15
n_s = n_cols*n_rows

n_sampled = 25 #number of sampled cells

N_surfaces = c(2,5,8)
Unbiased = c("all","one","none")  #whether generated surfaces are biased or not

n_sims = 100


XY = data.frame(lat = rep(c(1:n_cols),n_rows),lng= rep(c(1:n_rows),each=n_cols))
Dists = as.matrix(dist(XY, diag=T, upper=T))
Cov_exp = fields::Exp.cov(XY,XY,aRange=10)
L = chol(Cov_exp)
Cov_exp_bias = fields::Exp.cov(XY,XY,aRange=6)
L_bias = chol(Cov_exp_bias)

RMSE = array(NA,dim=c(length(N_surfaces),length(Unbiased),5,n_sims)) #next to last array dim is for estimator type

#set up SPDE basis using INLA
# Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
Grid_locs = expand.grid(y=c(1:n_rows),x=c(1:n_cols))
mesh = inla.mesh.create( Grid_locs )
n_knots = mesh$n
Eta_index = mesh$idx$loc-1  #which spde REs to apply as random effects for each cell centroid
#note it might be a lot faster to just to model the individual data points instead of the whole grid
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]

#initialize data with stuff that shouldn't change among simulation runs
Data1 <- list("M0"=spde$M0,"M1"=spde$M1,"M2"=spde$M2,"A_s"=rep(1,n_s),"P_i"=rep(1,n_sampled),
             "X_s"=matrix(1,n_s,1),"flag"=1,"n_i"=n_sampled,"n_s"=n_s,n_eta=n_knots,
             "Eta_index"=Eta_index,"options"=c(0,0),"beta_pri"=c(0,1),"matern_pri"=c(0,0,0,0))

Data_all <- list("M0"=spde$M0,"M1"=spde$M1,"M2"=spde$M2,
                 "flag"=1,"n_i"=n_sampled,"n_s"=n_s,n_mesh=n_knots,
                 "options"=c(0,0),"Mesh_index"=Eta_index,"matern_pri"=c(0,0,0,0))

Params1 <- list("Beta"=0, "log_tau"=0, "log_kappa"=0,"Etainput_s"=rep(0,n_knots))      
Random1 <- "Etainput_s"
Map1 <- NULL

Params_all <- list("log_tau_mu"=0, "log_kappa_mu"=0,"Mu_s"=rep(0,n_knots))      
Random_all <- c("Mu_s","Xi_s")

Ind_log = list("isim"=NA)
DEBUG=TRUE

st_time = Sys.time()
for(isim in 1:n_sims){
  Ind_log$isim=isim
  if(DEBUG)save(Ind_log,file="Ind_log.RData")
  set.seed(init.seed+isim)
  
  #simulate true spatial surface
  
  Mu = t(L) %*% rnorm(n_s)
  Nu = exp(Mu)
  Pi = Nu/(sum(Nu))
  
  #simulate biased spatial surfaces if applicable
  for(iscen in 1:length(N_surfaces)){
    Ind_log$iscen=iscen
    Mu_surf = Nu_surf = Pi_surf = matrix(0,n_s,N_surfaces[iscen])
    Data_all$n_surf=n_surf=N_surfaces[iscen]
    Params_all$Xi_s = matrix(0,n_knots,n_surf)
    Params_all$log_alpha = Params_all$log_tau_xi = Params_all$log_kappa_xi = rep(0,n_surf)

    for(ibias in 1:length(Unbiased)){
      Ind_log$ibias=ibias
      if(ibias==1){ #all surfaces biased
        for(isurf in 1:N_surfaces[iscen]){
          Mu_surf[,isurf] = Mu + t(L_bias) %*% rnorm(n_s)
          Nu_surf[,isurf] = exp(Mu_surf[,isurf])
          Pi_surf[,isurf] = Nu_surf[,isurf]/sum(Nu_surf[,isurf])
        }
      }
      if(ibias==2){ #one surface unbiased (we'll make it the first one)
        Mu_surf[,1]=Mu
        Pi_surf[,1] = exp(Mu)/sum(exp(Mu))
        for(isurf in 2:N_surfaces[iscen]){
          Mu_surf[,isurf] = Mu + t(L_bias) %*% rnorm(n_s)
          Nu_surf[,isurf] = exp(Mu_surf[,isurf])
          Pi_surf[,isurf] = Nu_surf[,isurf]/sum(Nu_surf[,isurf])
        }
      }
      if(ibias==3){ #no surface biased
        for(isurf in 1:N_surfaces[iscen]){
          Mu_surf[,isurf] = Mu 
          Nu_surf[,isurf] = Nu
          Pi_surf[,isurf] = Pi
        }
      }
      
      # simulate data (count?)
      Sampled = Count = matrix(0,n_sampled,N_surfaces[iscen])
      for(isurf in 1:N_surfaces[iscen]){
        Sampled[,isurf]=sample(c(1:n_s),n_sampled)
        Count[,isurf] = rpois(n_sampled,(N*Pi_surf[Sampled[,isurf],isurf]))
      }
      
      if(DEBUG)save(Ind_log,file="Ind_log.RData")
      
      # estimate surfaces for individual datasets; produce estimates of Y=log(Z) w/ Var-Cov matrix
      Y_ind = VC_Y_ind = vector("list",N_surfaces[iscen])
      for(isurf in 1:N_surfaces[iscen]){
        Ind_log$ind_surf=1
        Data1$C_i = Count[,isurf]
        Data1$S_i = Sampled[,isurf]-1
        
          # Make object
        #compile( paste0(Version,".cpp") )
        Start_time = Sys.time()
        #setwd( "C:/Users/paul.conn/git/OkhotskST/OkhotskSeal/src/")
        
        #st_time = Sys.time()
        
        Data1$options[1]=1  #use normalization trick
        Obj = MakeADFun( data=Data1, parameters=Params1, random=Random1, map=Map1, DLL="fit_1_surface",silent=FALSE)
        Obj$fn( Obj$par )
        Obj <- normalize ( Obj , flag ="flag", value = 0)
        Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
        Upper = 50
        Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
        Report = Obj$report()
        
        #cat(paste0("time = ",Sys.time()-st_time))
        
        SD=sdreport(Obj,bias.correct=FALSE)
        Y_ind[[isurf]] = SD$value
        VC_Y_ind[[isurf]] = SD$cov
        
        if(DEBUG)save(Ind_log,file="Ind_log.RData")
        
      }

      Ind_log$ind_surf=0
      Ind_log$inv=1
      det_flag = 0
      Data_all$Y_i = matrix(unlist(Y_ind),n_s,n_surf)
      Data_all$Omega_Y = array(dim=c(n_s,n_s,n_surf))   #inverse of Sigma needed for fast dmvnorm-type calcs
      for(isurf in 1:n_surf){
        Data_all$Omega_Y[,,isurf]=solve(VC_Y_ind[[isurf]])
        if(det(VC_Y_ind[[isurf]])<0 | is.na(det(VC_Y_ind[[isurf]])))det_flag=1
      }
      if(DEBUG)save(Ind_log,file="Ind_log.RData")
      
      Ind_log$inv=0
      
      # estimate combined surface assuming all surfaces biased
      if(det_flag==0){  #only run if all surfaces have pos def VC matrix
        Ind_log$est_mod=1
        if(DEBUG)save(Ind_log,file="Ind_log.RData")
        
        Data_all$options[1]=0  #use normalization trick
        Map_all=NULL
        Obj = MakeADFun( data=Data_all, parameters=Params_all, random=Random_all, map=Map_all, DLL="fit_multiple_surfaces",silent=FALSE)
        Obj$fn( Obj$par )
        #Obj <- normalize ( Obj , flag ="flag", value = 0)
        Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
        Upper = 50
        Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
        Report = Obj$report()
        #SD=sdreport(Obj,bias.correct=FALSE)
        
        RMSE[iscen,ibias,1,isim] = sqrt(sum((Pi-Report$Pi[1:n_s])^2)/n_s)
        
        
        # estimate combined surface assuming all surfaces unbiased
        Ind_log$est_mod=2
        if(DEBUG)save(Ind_log,file="Ind_log.RData")
        
        Map_all = list("Xi_s"=factor(matrix(NA,n_knots,n_surf)),"log_kappa_xi"=factor(rep(NA,n_surf)),
                       "log_tau_xi"=factor(rep(NA,n_surf)))
        Obj = MakeADFun( data=Data_all, parameters=Params_all, random=Random_all, map=Map_all, DLL="fit_multiple_surfaces",silent=FALSE)
        Obj$fn( Obj$par )
        #Obj <- normalize ( Obj , flag ="flag", value = 0)
        Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
        Upper = 50
        Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
        Report = Obj$report()
        #SD=sdreport(Obj,bias.correct=FALSE)
        
        RMSE[iscen,ibias,2,isim] = sqrt(sum((Pi-Report$Pi[1:n_s])^2)/n_s)
        
        
        # estimate combined surface assuming one surface unbiased
        Ind_log$est_mod=3
        if(DEBUG){
          save(Ind_log,file="Ind_log.RData")
          Inputs = list("Data"=Data_all, "Params"=Params_all,"Random"=Random_all,"Map_all"=Map_all)
          save(Inputs,file="Inputs_debug.RData")
        }
        
        Map_all = list("Xi_s"=factor(c(rep(NA,n_knots),c((2*n_surf-1):(n_knots*(n_surf-1)+2*n_surf-2)))),"log_kappa_xi"=factor(c(NA,c(1:(n_surf-1)))),
                       "log_tau_xi"=factor(c(NA,c(n_surf:(2*n_surf-2)))))
        cat("before Obj")
        Obj = MakeADFun( data=Data_all, parameters=Params_all, random=Random_all, map=Map_all, DLL="fit_multiple_surfaces",silent=FALSE)
        cat("after Obj")
        Obj$fn( Obj$par )
        #Obj <- normalize ( Obj , flag ="flag", value = 0)
        Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
        Upper = 50
        Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
        Report = Obj$report()
        #SD=sdreport(Obj,bias.correct=FALSE)
        
        RMSE[iscen,ibias,3,isim] = sqrt(sum((Pi-Report$Pi[1:n_s])^2)/n_s)
      }
      
      # just compute mean
      Pi_obs = matrix(0,n_surf,n_s)
      for(isurf in 1:n_surf)Pi_obs[isurf,]=exp(Data_all$Y_i[,isurf])/sum(exp(Data_all$Y_i[,isurf]))
      Pi_mean=colMeans(Pi_obs)
      RMSE[iscen,ibias,4,isim] = sqrt(sum((Pi-Pi_mean)^2)/n_s)
      
      
      # just take first surface (if it is the only unbiased one, is there any utility at including other surfaces?)
      RMSE[iscen,ibias,5,isim] = sqrt(sum((exp(Data_all$Y_i[,1])/sum(exp(Data_all$Y_i[,1]))-Pi)^2)/n_s)
      
    }
  }
  save(RMSE,file="sim_RMSE.RData")
  cat(paste0("Sim ",isim," Time elapsed: ", Sys.time()-st_time,"\n"))
}


#data mean
Pi_obs = matrix(0,n_surf,n_s)
for(isurf in 1:n_surf)Pi_obs[isurf,]=exp(Data_all$Y_i[,isurf])/sum(exp(Data_all$Y_i[,isurf]))
Pi_mean=colMeans(Pi_obs)

source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')

plot_N_map_xy(N=Pi,XY=Grid_locs,leg.title="Rel abundance")

plot_N_map_xy(N=Report$Pi[1:n_s],XY=Grid_locs,leg.title="Rel abundance")

plot_N_map_xy(N=Pi_mean,XY=Grid_locs,leg.title="Rel abundance")

sum((Pi-Report$Pi[1:n_s])^2)/n_s
sum((Pi-Pi_mean)^2)/n_s


#plot_N_map_xy(N=Report$Mu_s[1:n_s]/sum(Report$Mu_s[1:n_s]),XY=Grid_locs,leg.title="Rel abundance")

plot_N_map_xy(N=exp(Data_all$Y_i[,2]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Data_all$Y_i[,1]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Data_all$Y_i[,3]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Data_all$Y_i[,4]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Data_all$Y_i[,5]),XY=Grid_locs,leg.title="Rel abundance")

plot_N_map_xy(N=exp(Report$Xi[1:n_s,1]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Report$Xi[1:n_s,2]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Report$Xi[1:n_s,3]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Report$Xi[1:n_s,4]),XY=Grid_locs,leg.title="Rel abundance")
plot_N_map_xy(N=exp(Report$Xi[1:n_s,5]),XY=Grid_locs,leg.title="Rel abundance")




#cur_day=10
#plot_N_map_xy(N=Report$Z_s[1,((cur_day-1)*n_s+1):(cur_day*n_s)],XY=loc_s,leg.title="Abundance")

Converge=Opt$convergence


#plot count data by apparent species
Count = rowSums(Data$C_i[,1:3])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")
plot_N_map_xy(N=Report$E_count_sp[1,],XY=Coords3,leg.title="Count")

#look at thinning parameters
plot_N_map_xy(N=Thin[1:1285],XY=Coords3,leg.title="Count")
plot_N_map_xy(N=Report$Thin_trans[1,],XY=Coords3,leg.title="Count")

Count = rowSums(Data$C_i[,4:6])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

Count = rowSums(Data$C_i[,7:9])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

Count = rowSums(Data$C_i[,10:12])
crap = Count/Data$P_i
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=crap,XY=Coords3,leg.title="Count")
# 

Count = rowSums(Report$E_count_obs[,4:6])
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")
  
Count = Report$E_count_sp[2,]
Coords3 = Coords2[Mapping[,1],]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")

# 
#plot thinning by day
plot(DayHour[,1],Report$Thin_trans[2,])


#look at proportion of zeroes
get_pred_zeros <- function(E_count,sim=TRUE){
  if(sim)Rand = matrix(rpois(length(E_count),E_count),nrow(E_count),ncol(E_count))
  if(!sim)Rand = E_count
  Prop_zero=rep(0,ncol(E_count))
  for(i in 1:ncol(E_count))Prop_zero[i]=sum(Rand[,i]==0)
  Prop_zero
}

E_count = cbind(rowSums(Report$E_count_obs[,1:3]),rowSums(Report$E_count_obs[,4:6]),rowSums(Report$E_count_obs[,7:9]),rowSums(Report$E_count_obs[,10:12]))
Data_sp = cbind(rowSums(Data$C_i[,1:3]),rowSums(Data$C_i[,4:6]),rowSums(Data$C_i[,7:9]),rowSums(Data$C_i[,10:12]))
get_pred_zeros(E_count=E_count)
get_pred_zeros(Data_sp,sim=FALSE)


