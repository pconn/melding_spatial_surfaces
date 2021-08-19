# plot example simulation

#library( RandomFields )
library( TMB )
library( INLA )
#library( TMBhelper)  #install w/ devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
library( mvtnorm)
library(LaplacesDemon)
#library(TMBdebug)

jay_dir = 'c:/users/paul.conn/git/Data_integration/JaySimulationCode/'
source(paste0(jay_dir,'pointSimSyst.R'))
source(paste0(jay_dir,'pointSimCSR.R'))
source(paste0(jay_dir,'geostatSim.R'))
source(paste0(jay_dir,'distGeoAni.R'))
source(paste0(jay_dir,'corModels.R'))
source(paste0(jay_dir,'splmm.R'))
source(paste0(jay_dir,'covParmIni.R'))
source(paste0(jay_dir,'m2LLi.R'))
source(paste0(jay_dir,'m2LL.R'))
source(paste0(jay_dir,'makeCovMat.R'))
source(paste0(jay_dir,'addBreakColorLegend.R'))
source(paste0(jay_dir,'predict.splmm.R'))

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

#TmbFile_all = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_all_surfaces"
#compile(paste0(TmbFile_all,".cpp"),"-O1 -g",DLLFLAGS="") 

TmbFile1 = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_1_surface"
compile(paste0(TmbFile1,".cpp"),"-O1 -g",DLLFLAGS="") 

TmbFile2 = "C:/Users/paul.conn/git/Data_integration/Data_integration/src/fit_multiple_surfaces"
compile(paste0(TmbFile2,".cpp"),"-O1 -g",DLLFLAGS="") 

dyn.load( dynlib(TmbFile1) )
dyn.load( dynlib(TmbFile2) )


init.seed = 123456
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

isim=1
st_time = Sys.time()
  Ind_log$isim=isim
  if(DEBUG)save(Ind_log,file="Ind_log.RData")
  set.seed(init.seed+isim)
  
  #simulate true spatial surface
  
  Mu = t(L) %*% rnorm(n_s)
  Nu = exp(Mu)
  Pi = Nu/(sum(Nu))
  
  #simulate biased spatial surfaces if applicable
  
  iscen=2
    Ind_log$iscen=iscen
    Mu_surf = Nu_surf = Pi_surf = Xi_surf = matrix(0,n_s,N_surfaces[iscen])
    Data_all$n_surf=n_surf=N_surfaces[iscen]
    Params_all$Xi_s = matrix(0,n_knots,n_surf)
    Params_all$log_alpha = Params_all$log_tau_xi = Params_all$log_kappa_xi = rep(0,n_surf)
 
     #one surface unbiased (we'll make it the first one)
        Mu_surf[,1]=Mu
        Pi_surf[,1] = exp(Mu)/sum(exp(Mu))
        for(isurf in 2:N_surfaces[iscen]){
          Mu_surf[,isurf] = Mu + t(L_bias) %*% rnorm(n_s)
          Nu_surf[,isurf] = exp(Mu_surf[,isurf])
          Pi_surf[,isurf] = Nu_surf[,isurf]/sum(Nu_surf[,isurf])
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
        
      }
      
      # just compute mean
      Pi_obs = matrix(0,n_surf,n_s)
      for(isurf in 1:n_surf)Pi_obs[isurf,]=exp(Data_all$Y_i[,isurf])/sum(exp(Data_all$Y_i[,isurf]))
      Pi_mean=colMeans(Pi_obs)
 
      
      # just take first surface (if it is the only unbiased one, is there any utility at including other surfaces?)
 
    


#data mean
Pi_obs = matrix(0,n_surf,n_s)
for(isurf in 1:n_surf)Pi_obs[isurf,]=exp(Data_all$Y_i[,isurf])/sum(exp(Data_all$Y_i[,isurf]))
Pi_mean=colMeans(Pi_obs)

source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')

plot1 = plot_N_map_xy(N=N*Pi,XY=Grid_locs,leg.title="N")+
  scale_fill_viridis_c(limits = c(0,320), breaks = c(0, 100, 200, 300),name=expression(N(s)),values=c(0,.2,.4,1))
jpeg('True_abundance.jpg')
  plot1
dev.off()

plot2 = plot_N_map_xy(N=N*Pi_surf[,1],XY=Grid_locs,leg.title="N")+scale_fill_viridis_c(limits = c(0,320), breaks = c(0, 100, 200, 300),name="Biased N",values=c(0,.2,.4,1))

jpeg('true_bias1.jpg')
plot2
dev.off()

Plot_N = N*Pi_surf[,2]
Plot_N[which(Plot_N>250)]=250
plot3 = plot_N_map_xy(N=Plot_N,XY=Grid_locs,leg.title="Biased N")+scale_fill_viridis_c(limits = c(0,320), breaks = c(0, 100, 200, 300),name=expression("N"[i]*"(s)"),values=c(0,.2,.4,1))
jpeg('true_bias2.jpg')
plot3
dev.off()

plot4 = plot_N_map_xy(N=N*Pi_surf[,3],XY=Grid_locs,leg.title="Biased N")+scale_fill_viridis_c(limits = c(0,320), breaks = c(0, 100, 200, 300),name="Biased N",values=c(0,.2,.4,1))
jpeg('true_bias3.jpg')
plot4
dev.off()

plot5 = plot_N_map_xy(N=N*Pi_surf[,4],XY=Grid_locs,leg.title="Biased N")+scale_fill_viridis_c(limits = c(0,320), breaks = c(0, 100, 200, 300),name="Biased N",values=c(0,.2,.4,1))
jpeg('true_bias4.jpg')
plot5
dev.off()

Count_plot = rep(NA,n_s)
Count_plot[Sampled[,1]]=Count[,1]
plot6 = plot_N_map_xy(N=Count_plot,XY=Grid_locs,leg.title="Count")
jpeg('count1.jpg')
plot6
dev.off()

Count_plot = rep(NA,n_s)
Count_plot[Sampled[,2]]=Count[,2]
plot7 = plot_N_map_xy(N=Count_plot,XY=Grid_locs,leg.title=expression("C"[i]*"(s)"))
jpeg('count2.jpg')
plot7
dev.off()

Count_plot = rep(NA,n_s)
Count_plot[Sampled[,3]]=Count[,4]
plot8 = plot_N_map_xy(N=Count_plot,XY=Grid_locs,leg.title="Count")
jpeg('count3.jpg')
plot8
dev.off()

Count_plot = rep(NA,n_s)
Count_plot[Sampled[,4]]=Count[,4]
plot9 = plot_N_map_xy(N=Count_plot,XY=Grid_locs,leg.title="Count")
jpeg('count4.jpg')
plot9
dev.off()

Pi_est = Report$Pi[1:n_s]/sum(Report$Pi[1:n_s])
plot10 = plot_N_map_xy(N=Pi_est,XY=Grid_locs,leg.title="N")+
  scale_fill_viridis_c(limits = c(0,0.032), breaks = c(0, 0.01, 0.02, 0.03),name=expression(widehat(pi)*"(s)"),values=c(0,.2,.4,1))
jpeg('Combined_surface.jpg')
plot10
dev.off()


plot11 = plot_N_map_xy(N=exp(Report$Xi_s[1:n_s,1]),XY=Grid_locs,leg.title="Bias")+
  scale_fill_gradient2(low = scales::muted("red"),
                       mid = "white",
                       high = scales::muted("blue"),
                       midpoint = 1,name="Bias") 
jpeg('Xi1.jpg')
  plot11
dev.off()

plot12 = plot_N_map_xy(N=exp(Report$Xi_s[1:n_s,2]),XY=Grid_locs)+
  scale_fill_gradient2(low = scales::muted("blue"),
                       mid = "white",
                       high = scales::muted("red"),
                       midpoint = 1,name=expression(widehat(xi)[2]*"(s)"))
jpeg('Xi2.jpg')
plot12
dev.off()

plot13 = plot_N_map_xy(N=exp(Report$Xi_s[1:n_s,3]),XY=Grid_locs,leg.title="Bias")+
  scale_fill_gradient2(low = scales::muted("blue"),
                       mid = "white",
                       high = scales::muted("red"),
                       midpoint = 1,name=expression(widehat(xi)[3]*"(s)")) 
jpeg('Xi3.jpg')
plot13
dev.off()

plot14 = plot_N_map_xy(N=exp(Report$Xi_s[1:n_s,4]),XY=Grid_locs,leg.title="Bias")+
  scale_fill_gradient2(low = scales::muted("blue"),
                       mid = "white",
                       high = scales::muted("red"),
                       midpoint = 1,name=expression(widehat(xi)[4]*"(s)")) 
jpeg('Xi4.jpg')
plot14
dev.off()

plot15 = plot_N_map_xy(N=exp(Data_all$Y_i[,2]),XY=Grid_locs,leg.title="Est. N")+
  scale_fill_viridis_c(limits = c(0,320), breaks = c(0, 100, 200, 300),name=expression(widehat(N)*"(s)"),values=c(0,.2,.4,1))

plot16 = plot_N_map_xy(N=exp(Data_all$Y_i[,1])/sum(exp(Data_all$Y_i[,1])),XY=Grid_locs)+
  scale_fill_viridis_c(limits = c(0,0.032), breaks = c(0, 0.01, 0.02, 0.03),name=expression(widehat(pi)[1]*"(s)"),values=c(0,.2,.4,1))

Cur_pi = exp(Data_all$Y_i[,2])/sum(exp(Data_all$Y_i[,2]))
Cur_pi[which(Cur_pi>0.025)]=0.025
plot17 = plot_N_map_xy(N=Cur_pi,XY=Grid_locs)+
  scale_fill_viridis_c(limits = c(0,0.032), breaks = c(0, 0.01, 0.02, 0.03),name=expression(widehat(pi)[2]*"(s)"),values=c(0,.2,.4,1))

plot18 = plot_N_map_xy(N=exp(Data_all$Y_i[,3])/sum(exp(Data_all$Y_i[,3])),XY=Grid_locs)+
  scale_fill_viridis_c(limits = c(0,0.032), breaks = c(0, 0.01, 0.02, 0.03),name=expression(widehat(pi)[3]*"(s)"),values=c(0,.2,.4,1))

plot19 = plot_N_map_xy(N=exp(Data_all$Y_i[,4])/sum(exp(Data_all$Y_i[,4])),XY=Grid_locs)+
  scale_fill_viridis_c(limits = c(0,0.032), breaks = c(0, 0.01, 0.02, 0.03),name=expression(widehat(pi)[4]*"(s)"),values=c(0,.2,.4,1))


plot1 = plot1 + theme(plot.title = element_text(hjust = 0))+ggtitle('A.')
plot3 = plot3 + theme(plot.title = element_text(hjust = 0))+ggtitle('B.')
plot7 = plot7 + theme(plot.title = element_text(hjust = 0))+ggtitle('C.')
plot15 = plot15 + theme(plot.title = element_text(hjust = 0))+ggtitle('D.')

pdf("Sim_plot1.pdf",height=6,width=8)
gridExtra::grid.arrange(plot1,plot3,plot7,plot15,ncol=2)
#cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()

plot16 = plot16 + theme(plot.title = element_text(hjust = 0))+ggtitle('A.') 
plot10 = plot10 + theme(plot.title = element_text(hjust = 0))+ggtitle('B.')
plot17 = plot17 + theme(plot.title = element_text(hjust = 0))+ggtitle('C.')
plot12 = plot12 + theme(plot.title = element_text(hjust = 0))+ggtitle('D.')
plot18 = plot18 + theme(plot.title = element_text(hjust = 0))+ggtitle('E.')
plot13 = plot13 + theme(plot.title = element_text(hjust = 0))+ggtitle('F.')
plot19 = plot19 + theme(plot.title = element_text(hjust = 0))+ggtitle('G.')
plot14 = plot14 + theme(plot.title = element_text(hjust = 0))+ggtitle('H.')

###plot 2: surfaces ---> estimates
pdf("Sim_plot2.pdf",height=6,width=8)
gridExtra::grid.arrange(plot16,plot10,plot17,plot12,plot18,plot13,plot19,plot14,ncol=2)
#cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()




