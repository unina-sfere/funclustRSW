# fclust ------------------------------------------------------------------

fit_fclust<-function(data,...){
  
  
  n_obs<-dim(data$X)[2]
  data_vec<-list(x=matrix(data$X,ncol=1),timeindex=rep(1:length(data$grid),n_obs),curve=rep(1:n_obs,each=length(data$grid)))
  
  testfit <- fitfclust(data=data_vec,...)
}



fit_fclust_ms<-function(data,num_cluster_seq=2:4,dim_seq=5,ncores=1){
  
  comb_list<-expand.grid(dim_seq,num_cluster_seq)
  comb_list<-cbind(comb_list,comb_list[,2]-1)
  n_obs <- ncol(data$X)
  data_vec<-list(x = matrix(data$X,ncol=1),
                 timeindex = rep(1:length(data$grid),n_obs),
                 curve = rep(1:n_obs,each=length(data$grid)))
  N<-length(unique(data_vec$curve))
  
  parr_fun<-function(ii){
    cat(ii)
    parameters<-as.numeric(comb_list[ii,])
    dim_i<-parameters[1]
    num_cluster<-parameters[2]
    h<-parameters[3]
    fit_fclust<-fitfclust(data=data_vec,trace = TRUE,K=num_cluster,q  =dim_i,h=h,pert =1,grid = data$grid ,plot=TRUE)
    ll<-loglik(data=fit_fclust$data, parameters=fit_fclust$parameters, vars=fit_fclust$vars, FullS=fit_fclust$FullS)
    N<-length(unique(data_vec$curve))
    K <- dim(fit_fclust$vars$gamma)[2]
    q <- dim(fit_fclust$vars$gamma)[3]
    h<- dim(fit_fclust$parameters$alpha)[2]
    df<-(K-1)+1+(q^2+q)/2+q+K*h+h*q-h-h^2
    BIC<- -2*ll+log(N)*df
    AIC<- -2*ll+2*df
    piigivej <- fit_fclust$vars$piigivej
    class <- apply(piigivej, 1, which.max)
    
    parameters <- fit_fclust$parameters
    FullS <- fit_fclust$FullS
    cluster.mean <- matrix(0,K,dim(FullS)[1])
    for(k in 1:K)
      cluster.mean[k,] <- FullS %*% (parameters$lambda.zero + parameters$Lambda %*% parameters$alpha[k,  ])
    
    out<-list(AIC=AIC,
              BIC=BIC,
              piigivej=piigivej,
              class=class,
              cluster.mean=cluster.mean)
  }
  
  
  
  cores <- detectCores()
  vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_fun,mc.cores = ncores)
  AIC<-sapply(vec_par,"[[",1)
  BIC<-sapply(vec_par,"[[",2)
  piigivej<-lapply(vec_par,"[[",3)
  class<-lapply(vec_par,"[[",4)
  cluster.mean<-lapply(vec_par,"[[",5)
  par_opt_AIC<-as.numeric(comb_list[max(which(AIC==min(AIC))),])
  par_opt_BIC<-as.numeric(comb_list[max(which(BIC==min(BIC))),])
  piigivej_opt_AIC <- piigivej[[max(which(AIC==min(AIC)))]]
  piigivej_opt_BIC <- piigivej[[max(which(BIC==min(BIC)))]]
  class_opt_AIC <- class[[max(which(AIC==min(AIC)))]]
  class_opt_BIC <- class[[max(which(BIC==min(BIC)))]]
  cluster.mean_opt_AIC <- cluster.mean[[max(which(AIC==min(AIC)))]]
  cluster.mean_opt_BIC <- cluster.mean[[max(which(BIC==min(BIC)))]]
  
  
  out<-list(par_opt_AIC=par_opt_AIC,
            par_opt_BIC=par_opt_BIC,
            AIC=AIC,
            BIC=BIC,
            comb_list=comb_list,
            piigivej = piigivej,
            class = class,
            par_opt_AIC = par_opt_AIC,
            par_opt_BIC = par_opt_BIC,
            piigivej_opt_AIC = piigivej_opt_AIC,
            piigivej_opt_BIC = piigivej_opt_BIC,
            class_opt_AIC = class_opt_AIC,
            class_opt_BIC = class_opt_BIC,
            cluster.mean=cluster.mean,
            cluster.mean_opt_AIC=cluster.mean_opt_AIC,
            cluster.mean_opt_BIC=cluster.mean_opt_BIC
            )
  return(out)
}



# Curvclust ---------------------------------------------------------------
curvclust<-function(data,K, structure, mixed, reduction){
  
  N<-ncol(data$X)
  q<-nrow(data$X)
  fdat= list()
  for (j in 1:N) fdat[[j]] =data$X[,j]
  CCD = new("CClustData",Y=fdat,filter.number=1)
  
  if (reduction) CCDred = getUnionCoef(CCD) else CCDred <- CCD
  
  CCO = new("CClustO")
  CCO@nbclust= K
  if (K > 1) {
    CCO@Gamma2.structure = structure
    CCO@loglikplot=TRUE
    CCO@init="SEM"
    CCO@burn=100
  }
  
  CCR <- if (K == 1) {
    getFMM(CCDred, CCO) 
    } else {
      if (mixed) getFCMM(CCDred,CCO) else getFCM(CCDred,CCO)
      }
  
  mu<-matrix(unlist(getwr.mu(CCR,CCO,CCDred)),q,K)
  
  matplot(mu,type="l")
  cluster <- if (K > 1) apply(CCR@Tau,1,which.max) else NULL
  
  BIC<--getBIC(CCR,CCDred)
  ICL<-getICL(CCR,CCDred)
  
  out<-list(CCR=CCR,
            mu=mu,
            BIC=BIC,
            ICL=ICL,
            cluster = cluster)
  return(out)
}
curvclust_ms<-function(data,num_cluster_seq=2:3,ncores=1,
                       structure,mixed, reduction){
  
  if(mixed == FALSE & structure != "none") 
    stop(print("if mixed FALSE only none is possible"))
  
  if(mixed == TRUE & structure == "none") 
    stop(print("none not possible with mixed"))
  
  parr_fun<-function(ii){
    cat(ii)
    
    num_cluster<-num_cluster_seq[ii]
    
    fit_curvclust<-curvclust(data=data,K=num_cluster,structure=structure,
                             mixed = mixed, reduction = reduction)
    
    BIC<- fit_curvclust$BIC
    # ICL<- fit_curvclust$ICL
    mu <- fit_curvclust$mu
    tau <- fit_curvclust$CCR@Tau
    class <- if (num_cluster > 1) apply(fit_curvclust$CCR@Tau, 1, which.max) else
      rep(1, ncol(data$X))
    
    out<-list(
      # ICL=ICL,
              BIC=BIC,
              tau=tau,
              class=class,
              mu=mu)
  }
  
  cores <- detectCores()
  vec_par<-mclapply(seq(1,length(num_cluster_seq)),parr_fun,mc.cores = ncores)
  # ICL<-sapply(vec_par,"[[",1)
  BIC<-sapply(vec_par,"[[",1)
  tau<-lapply(vec_par,"[[",2)
  class<-lapply(vec_par,"[[",3)
  mu<-lapply(vec_par,"[[",4)
  # par_opt_ICL<-num_cluster_seq[max(which(ICL==min(ICL, na.rm = TRUE)))]
  par_opt_BIC<-num_cluster_seq[max(which(BIC==max(BIC, na.rm = TRUE)))]
  # tau_opt_ICL <- tau[[max(which(ICL==min(ICL, na.rm = TRUE)))]]
  tau_opt_BIC <- tau[[max(which(BIC==max(BIC, na.rm = TRUE)))]]
  # class_opt_ICL <- class[[max(which(ICL==min(ICL, na.rm = TRUE)))]]
  class_opt_BIC <- class[[max(which(BIC==max(BIC, na.rm = TRUE)))]]
  # mu_opt_ICL <- mu[[max(which(ICL==min(ICL, na.rm = TRUE)))]]
  mu_opt_BIC <- mu[[max(which(BIC==max(BIC, na.rm = TRUE)))]]
  
  
  out<-list(
    # par_opt_ICL=par_opt_ICL,
            par_opt_BIC=par_opt_BIC,
            # ICL=ICL,
            BIC=BIC,
            tau=tau,
            class=class,
            mu=mu,
            tau_opt_BIC=tau_opt_BIC,
            # tau_opt_ICL=tau_opt_ICL,
            class_opt_BIC=class_opt_BIC,
            # class_opt_ICL=class_opt_ICL,
            # mu_opt_ICL=mu_opt_ICL,
            mu_opt_BIC=mu_opt_BIC)
}





# funHDDC -----------------------------------------------------------------

slopeHeuristic_mod<-function (mod) 
{
  main_data = mod$allCriteria
  who_notNA_norInfinite = !is.na(main_data$complexity) & is.finite(main_data$LL)
  n_valid = sum(who_notNA_norInfinite)
  if (n_valid == 0) {
    stop("There is not any valid model to be selected.")
  }
  else if (n_valid <= 2) {
    stop("At least 3 valid models are necessary to perform the slope heuristic. Otherwise, use another criterion.")
  }
  main_data = main_data[who_notNA_norInfinite, ]
  fit = MASS::rlm(LL ~ jitter(complexity), data = main_data, method = "MM")
  fit_coef = fit$coefficients
  if (fit_coef[2] < 0) {
    fit_coef[2] = 0
  }
  llpen = main_data$LL - 2 * fit_coef[2] * main_data$complexity
  SH = 2 * llpen
  
  plot(main_data$complexity, main_data$LL, type = "p", 
       xlab = "Model dimension", ylab = "Log-likelihood")
  abline(fit, col = "red")
  plot(main_data$rank, llpen, type = "p", xlab = "Rank of the model", 
       ylab = "Penalized log-likelihood")
  points(main_data$rank[which.max(llpen)], max(llpen), 
         pch = 19, col = "red")
  return(c(main_data$rank[which.max(llpen)],SH[which.max(llpen)]))
  
}


fit_funHDDC_ms<-function(data,num_cluster_seq=2:3,threshold_seq=seq(0.1,0.99,length.out = 3),
                         model = c('AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'),
                         ncores = 1,
                         nb.rep){
  
  basis<- create.bspline.basis(c(0,1), nbasis=20)
  loglam         = seq(-10, -4, 0.25)
  Gcvsave        = numeric()
  for(i in 1:length(loglam)){
    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    Sm.i    = smooth.basis(data$grid, data$X, fdPari)
    Gcvsave[i] = sum(Sm.i$gcv)
    
  }
  
  # Figure 5.1.
  lambda_s=10^loglam[which.min(Gcvsave)]
  
  plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
       ylab=expression(GCV(lambda)), lwd=2 )
  
  fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
  X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  mod<-funHDDC(X_fd,num_cluster_seq,
               model=model,
               show=TRUE,
               threshold=threshold_seq,
               kmeans.control = list(nstart=100, iter.max = 100),
               nb.rep = nb.rep)
  
  mu_fd<-fd(t(mod$mu),basis)
  
  SH<- slopeHeuristic_mod(mod)
  K_opt<-mod$allCriteria[SH[1],]$K 
  thr_opt<-mod$allCriteria[SH[1],]$threshold
  mod_opt<-funHDDC(X_fd,K_opt,model=c("AkjBkQkDk"),show=TRUE,
                   threshold=thr_opt,kmeans.control = list(nstart=100))
  out<-list(mod=mod,mod_opt=mod_opt,par_opt=c(lambda_s,K_opt,thr_opt,SH[2]))
  
  return(out)
}





# Distance-based ----------------------------------------------------------

distance_ss<-function(data,num_cluster=3){
  
  if(class(data)!="fd"){
    basis<- create.bspline.basis(c(0,1), nbasis=20)
    loglam         = seq(-10, -4, 0.25)
    Gcvsave        = numeric()
    for(i in 1:length(loglam)){
      fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
      Sm.i    = smooth.basis(data$grid, data$X, fdPari)
      Gcvsave[i] = sum(Sm.i$gcv)
      
    }
    
    
    lambda_s=10^loglam[which.min(Gcvsave)]
    
    plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
         ylab=expression(GCV(lambda)), lwd=2 )
    
    fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
    X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  }
  else{
    X_fd=data
  }
  X_fdata<-fdata(X_fd)
  mod<-kmeans.fd(X_fdata,ncl = num_cluster,par.metric = list(lp=2),par.dfunc = list(trim=0),
                 draw=FALSE, cluster.size = 1, max.iter = 5000)
  dist<-dist(metric.lp(X_fdata,p=2))
  clus<-mod$cluster
  mean_fd<-fdata2fd(mod$centers)
  out<-list(clus=clus,
            mean_fd=mean_fd,
            dist=dist)
  return(out)
  
}

distance_ms<-function(data,num_cluster_seq=2:7,met="sil", ncores = 1,
                      neighbSize=50){
  
  basis<- create.bspline.basis(c(0,1), nbasis=20)
  loglam         = seq(-10, -4, 0.25)
  Gcvsave        = numeric()
  for(i in 1:length(loglam)){
    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    Sm.i    = smooth.basis(data$grid, data$X, fdPari)
    Gcvsave[i] = sum(Sm.i$gcv)
    
  }
  
  
  lambda_s=10^loglam[which.min(Gcvsave)]
  
  plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
       ylab=expression(GCV(lambda)), lwd=2 )
  
  fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
  X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  
  parr_fun<-function(lll){
    print(lll)
    num_cluster_i<-num_cluster_seq[lll]
    
    mod<-distance_ss(X_fd,num_cluster_i)
    con<-connectivity(distance = mod$dist,mod$clus, neighbSize = neighbSize)
    dunn<-dunn(mod$dist,mod$clus)
    sil<-silhouette(mod$clus,mod$dist)
    # sil_median <- median(sil[,3])
    sil<-summary(sil)$avg.width
    out<-list(con=con,
              dunn=dunn,
              sil=sil)
              # sil_median=sil_median)
    return(out)
    
  }
  
  vec_par<-mclapply(seq(1,length(num_cluster_seq)),parr_fun, mc.cores = ncores)
  con<-sapply(vec_par,"[[",1)
  dunn<-sapply(vec_par,"[[",2)
  sil<-sapply(vec_par,"[[",3)
  # sil_median<-sapply(vec_par,"[[",4)
  
  con_opt<-as.numeric(num_cluster_seq[min(which(con==min(con)))])
  dunn_opt<-as.numeric(num_cluster_seq[min(which(dunn==max(dunn)))])
  sil_opt<-as.numeric(num_cluster_seq[min(which(sil==max(sil)))])
  # sil_median_opt<-as.numeric(num_cluster_seq[min(which(sil==max(sil_median)))])
  
  
  mod_opt_con<-distance_ss(X_fd,con_opt)
  mod_opt_dunn<-distance_ss(X_fd,dunn_opt)
  mod_opt_sil<-distance_ss(X_fd,sil_opt)
  # mod_opt_sil_medina<-distance_ss(X_fd,sil_opt_median)
  
  if(met=="sil"){
    out<-list(mod_opt_sil=mod_opt_sil,
              sil_opt=sil_opt,
              sil=sil)
  }
  else if(met=="con"){
    out<-list(mod_opt_con=mod_opt_con,
              con_opt=con_opt,
              con=con)
  }
  else if(met=="dunn"){
    out<-list(mod_opt_dunn=mod_opt_dunn,
              dunn_opt=dunn_opt,
              dunn=dunn)
  }
  else{
    out<-list(mod_opt_con=mod_opt_con,
              mod_opt_dunn=mod_opt_dunn,
              mod_opt_sil=mod_opt_sil,
              con_opt=con_opt,
              dunn_opt=dunn_opt,
              sil_opt=sil_opt,
              con=con,
              dunn=dunn,
              sil=sil)
  }
  return(out)
  
}

# Filtering B-spline ------------------------------------------------------


fil_bspline_ms2<-function(data,num_cluster_seq=2:5,grid=NULL){
  
  basis<- create.bspline.basis(c(0,1), nbasis=60)
  loglam         = seq(-10, -4, 0.25)
  Gcvsave        = numeric()
  for(i in 1:length(loglam)){
    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    Sm.i    = smooth.basis(data$grid, data$X, fdPari)
    Gcvsave[i] = sum(Sm.i$gcv)
    
  }
  
  
  lambda_s=10^loglam[which.min(Gcvsave)]
  
  plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
       ylab=expression(GCV(lambda)), lwd=2 )
  
  fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
  X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  
  B<-t(X_fd$coefs)
  gap_hc<-fviz_nbclust(B, hcut, method = "gap",k.max=num_cluster_seq[length(num_cluster_seq)])$data$gap
  gap_km<-fviz_nbclust(B, kmeans, method = "gap",k.max=num_cluster_seq[length(num_cluster_seq)])$data$gap
  model<-Mclust(B,G=num_cluster_seq)
  
  model$mean_fd <- fd(coef = model$parameters$mean, basis = X_fd$basis)
  
  K_opt_hc<-(1:num_cluster_seq[length(num_cluster_seq)])[which.max(gap_hc)]
  K_opt_km<-(1:num_cluster_seq[length(num_cluster_seq)])[which.max(gap_km)]
  
  mod_opt_hc<-hcut(B,k =K_opt_hc )
  mod_opt_km<-kmeans(B,centers =K_opt_km )
  
  mod_opt_mod<-model
  
  mean_fd<-list()
  mean_fd_hc<-fd(sapply(1:K_opt_hc, function(ll)colMeans(B[mod_opt_hc$cluster==ll,])),basis)
  mean_fd_km<-fd(t(mod_opt_km$centers),basis)
  mean_fd_model<-fd(mod_opt_mod$parameters$mean,basis)
  
  
  out<-list(mod_opt=list(ind_hc = mod_opt_hc,ind_km = mod_opt_km,mod_opt = mod_opt_mod),
            mean_fd=list(hc = mean_fd_hc,km = mean_fd_km,model = mean_fd_model),
            gap_hc=gap_hc,
            gap_km=gap_km)
  
  return(out)
  
}

fil_bspline_ms<-function(data,num_cluster_seq=2:5,grid=NULL,met="sil"){
  
  basis<- create.bspline.basis(c(0,1), nbasis=20)
  # loglam         = seq(-10, -4, 0.25)
  # Gcvsave        = numeric()
  # for(i in 1:length(loglam)){
  #   fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
  #   Sm.i    = smooth.basis(data$grid, data$X, fdPari)
  #   Gcvsave[i] = sum(Sm.i$gcv)
  #   
  # }
  # 
  # 
  # lambda_s=10^loglam[which.min(Gcvsave)]
  # 
  # plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
  #      ylab=expression(GCV(lambda)), lwd=2 )
  
  # fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
  
  fdPari  = fdPar(basis, Lfdobj=2,0)
  X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  
  B<-t(X_fd$coefs)
  mod<-clValid(B, num_cluster_seq, clMethods=c("hierarchical","kmeans"),
               validation="internal",method = "ward",neighbSize = 30,
               maxitems = ncol(data$X))
  model<-Mclust(B,G=c(1,num_cluster_seq))
  
  model$mean_fd <- fd(coef = model$parameters$mean, basis = X_fd$basis)
  
  con<-mod@measures[1,,]
  dunn<-mod@measures[2,,]
  sil<-mod@measures[3,,]
  K_opt_con<-apply(con, 2, which.min)
  K_opt_dunn<-apply(dunn, 2, which.max)
  K_opt_sil<-apply(sil, 2, which.max)
  K_opt<-list(K_opt_con,K_opt_dunn,K_opt_sil)
  mod_opt_hc<-lapply(1:3,function(ii)cutree(mod@clusterObjs$hierarchical,num_cluster_seq[K_opt[[ii]][1]]))
  mod_opt_km<-lapply(1:3,function(ii)mod@clusterObjs$kmeans[[K_opt[[ii]][2]]])
  mod_opt_mod<-model
  
  mean_fd<-list()
  mean_fd_hc<-lapply(1:3,function(ii)fd(sapply(1:length(unique(mod_opt_hc[[ii]])), function(ll)colMeans(B[mod_opt_hc[[ii]]==ll,])),basis))
  mean_fd_km<-lapply(1:3,function(ii)fd(t(mod_opt_km[[ii]]$centers),basis))
  mean_fd_model<-fd(mod_opt_mod$parameters$mean,basis)
  
  if(met=="sil"){
    ind=3
    out<-list(mod_opt_sil=list(ind_hc = mod_opt_hc[[ind]],ind_km = mod_opt_km[[ind]],mod_opt = mod_opt_mod),
              mean_fd=list(hc = mean_fd_hc[[ind]],km = mean_fd_km[[ind]],model = mean_fd_model),
              sil=sil)
  }
  else if(met=="con"){
    ind=1
    out<-list(mod_opt_con=list(ind_hc = mod_opt_hc[[ind]],ind_km = mod_opt_km[[ind]],mod_opt = mod_opt_mod),
              mean_fd=list(hc = mean_fd_hc[[ind]],km = mean_fd_km[[ind]],model = mean_fd_model),
              con=con)
  }
  else if(met=="dunn"){
    ind=2
    out<-list(mod_opt_dunn=list(ind_hc = mod_opt_hc[[ind]],ind_km = mod_opt_km[[ind]],mod_opt = mod_opt_mod),
              mean_fd=list(hc= mean_fd_hc[[ind]],km = mean_fd_km[[ind]],model = mean_fd_model),
              dunn=dunn)
  }
  else{
    out<-list(mod_opt_con=mod_opt_con,
              mod_opt_dunn=mod_opt_dunn,
              mod_opt_sil=mod_opt_sil,
              mean_fd_hc=mean_fd_hc,
              mean_fd_km=mean_fd_km,
              mean_fd_model=mean_fd_model,
              con=con,
              dunn=dunn,
              sil=sil)
  }
  return(out)
  
}
fil_bspline_ms_nclust<-function(data,num_cluster_seq=2:5,grid=NULL,nbasis=20){
  
  basis<- create.bspline.basis(c(0,1), nbasis=nbasis)
  loglam         = seq(-10, -4, 0.25)
  Gcvsave        = numeric()
  for(i in 1:length(loglam)){
    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    Sm.i    = smooth.basis(data$grid, data$X, fdPari)
    Gcvsave[i] = sum(Sm.i$gcv)

  }


  lambda_s=10^loglam[which.min(Gcvsave)]

  plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
       ylab=expression(GCV(lambda)), lwd=2 )
  
  # fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
  
  fdPari  = fdPar(basis, Lfdobj=2,0)
  X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  
  B<-t(X_fd$coefs)
  
  mod_km<-NbClust(data = B,  method = "kmeans",min.nc = 2, max.nc = num_cluster_seq[length(num_cluster_seq)])
  mod_hc<-NbClust(data = B,  method = "ward.D2",min.nc = 2, max.nc = num_cluster_seq[length(num_cluster_seq)])
  
  
  mod_km_kopt<-(mod_km$Best.nc[1,])[!mod_km$Best.nc[1,]==0]
  mod_hc_kopt<-(mod_hc$Best.nc[1,])[!mod_hc$Best.nc[1,]==0]
  par(mfrow=c(1,2))
  hist(mod_km_kopt,breaks = 0:num_cluster_seq[length(num_cluster_seq)])
  hist(mod_hc_kopt,breaks = 0:num_cluster_seq[length(num_cluster_seq)])
  K_opt_km<-as.numeric(names(sort(table(mod_km_kopt),decreasing=TRUE))[1])
  K_opt_hc<-as.numeric(names(sort(table(mod_hc_kopt),decreasing=TRUE))[1])
  
 
  ciao <- mclustBIC(B,G=c(1,num_cluster_seq),
                    initialization = list(subset = sample(1:ncol(data$X), as.integer(ncol(data$X) / 3))))
  plot(ciao)
  
  ciao <- mclustICL(B,G=c(1,num_cluster_seq))
  plot(ciao)
  
  model<-Mclust(B,G=c(1,num_cluster_seq), 
                initialization = list(subset = sample(1:ncol(data$X), as.integer(ncol(data$X) / 3))))

  mod_opt_hc<-hcut(B,k =K_opt_hc )
  mod_opt_km<-kmeans(B,centers =K_opt_km )
  
  mod_opt_mod<-model
  
  mean_fd<-list()
  mean_fd_hc<-fd(sapply(1:K_opt_hc, function(ll)colMeans(B[mod_opt_hc$cluster==ll,])),basis)
  mean_fd_km<-fd(t(mod_opt_km$centers),basis)
  mean_fd_model<-fd(mod_opt_mod$parameters$mean,basis)
  
  
  out<-list(mod_opt=list(ind_hc = mod_opt_hc,ind_km = mod_opt_km,mod_opt = mod_opt_mod),
            mean_fd=list(hc = mean_fd_hc,km = mean_fd_km,model = mean_fd_model),
            mod_km=mod_km,
            mod_hc=mod_hc)
  return(out)
  
}

# fil_bspline_ms<-function(data,num_cluster_seq=2:5,grid=NULL,met="sil"){
#   
#   
#   basis_vec         = seq(5,40,by=1)
#   Gcvsave        = numeric()
#   for(i in 1:length(basis_vec)){
#     
#     basis<- create.bspline.basis(c(0,1), nbasis=basis_vec[i])
#     Sm.i    = smooth.basis(data$grid, data$X, basis)
#     Gcvsave[i] = sum(Sm.i$gcv)
#     
#   }
#   
#   
#   basis_opt=basis_vec[which.min(Gcvsave)]
#   
#   plot(basis_vec, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
#        ylab=expression(GCV(lambda)), lwd=2 )
#   
#   basis<- create.bspline.basis(c(0,1), nbasis=basis_opt)
#   
#   X_fd<-smooth.basis(data$grid, data$X, basis)$fd
#   
#   B<-t(X_fd$coefs)
#   mod<-clValid(B, num_cluster_seq, clMethods=c("hierarchical","kmeans"),
#                validation="internal",method = "ward",neighbSize = 30,
#                maxitems = ncol(data$X))
#   model<-Mclust(B,G=num_cluster_seq)
#   
#   model$mean_fd <- fd(coef = model$parameters$mean, basis = X_fd$basis)
#   
#   con<-mod@measures[1,,]
#   dunn<-mod@measures[2,,]
#   sil<-mod@measures[3,,]
#   K_opt_con<-apply(con, 2, which.min)
#   K_opt_dunn<-apply(dunn, 2, which.max)
#   K_opt_sil<-apply(sil, 2, which.max)
#   K_opt<-list(K_opt_con,K_opt_dunn,K_opt_sil)
#   mod_opt_hc<-lapply(1:3,function(ii)cutree(mod@clusterObjs$hierarchical,num_cluster_seq[K_opt[[ii]][1]]))
#   mod_opt_km<-lapply(1:3,function(ii)mod@clusterObjs$kmeans[[K_opt[[ii]][2]]])
#   mod_opt_mod<-model
#   
#   mean_fd<-list()
#   mean_fd_hc<-lapply(1:3,function(ii)fd(sapply(1:length(unique(mod_opt_hc[[ii]])), function(ll)colMeans(B[mod_opt_hc[[ii]]==ll,])),basis))
#   mean_fd_km<-lapply(1:3,function(ii)fd(t(mod_opt_km[[ii]]$centers),basis))
#   mean_fd_model<-fd(mod_opt_mod$parameters$mean,basis)
#   
#   if(met=="sil"){
#     ind=3
#     out<-list(mod_opt_sil=list(ind_hc = mod_opt_hc[[ind]],ind_km = mod_opt_km[[ind]],mod_opt = mod_opt_mod),
#               mean_fd=list(hc = mean_fd_hc[[ind]],km = mean_fd_km[[ind]],model = mean_fd_model),
#               sil=sil)
#   }
#   else if(met=="con"){
#     ind=1
#     out<-list(mod_opt_con=list(ind_hc = mod_opt_hc[[ind]],ind_km = mod_opt_km[[ind]],mod_opt = mod_opt_mod),
#               mean_fd=list(hc = mean_fd_hc[[ind]],km = mean_fd_km[[ind]],model = mean_fd_model),
#               con=con)
#   }
#   else if(met=="dunn"){
#     ind=2
#     out<-list(mod_opt_dunn=list(ind_hc = mod_opt_hc[[ind]],ind_km = mod_opt_km[[ind]],mod_opt = mod_opt_mod),
#               mean_fd=list(hc= mean_fd_hc[[ind]],km = mean_fd_km[[ind]],model = mean_fd_model),
#               dunn=dunn)
#   }
#   else{
#     out<-list(mod_opt_con=mod_opt_con,
#               mod_opt_dunn=mod_opt_dunn,
#               mod_opt_sil=mod_opt_sil,
#               mean_fd_hc=mean_fd_hc,
#               mean_fd_km=mean_fd_km,
#               mean_fd_model=mean_fd_model,
#               con=con,
#               dunn=dunn,
#               sil=sil)
#   }
#   return(out)
#   
# }



# Filtering FPCA ----------------------------------------------------------



fil_fpca_ss_nbclust<-function(data,num_cluster_seq=3:4,per_comp=0.9,grid=NULL,...){
  
  
  basis<- create.bspline.basis(c(0,1), nbasis=100)
  loglam         = seq(-10, -4, 0.25)
  Gcvsave        = numeric()
  for(i in 1:length(loglam)){
    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    Sm.i    = smooth.basis(data$grid, data$X, fdPari)
    Gcvsave[i] = sum(Sm.i$gcv)
    
  }
  
  
  lambda_s=10^loglam[which.min(Gcvsave)]
  
  plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
       ylab=expression(GCV(lambda)), lwd=2 )
  
  fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
  X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  
  
  # sd_fd<-sd.fd(X_fd)
  mean_fd_ini<-mean.fd(X_fd)
  # X_std<-center.fd(X_fd)*(sd_fd)^-1
  princomp<-pca.fd(X_fd,min(X_fd$basis$nbasis, 30))
  # cum<-cumsum(princomp$values)/(X_fd$basis$rangeval[2]-X_fd$basis$rangeval[1])
  cum<-cumsum(princomp$values)/sum(princomp$values)
  num_com<-min(which(cum>per_comp))
  princomp<-pca.fd(X_fd,num_com)
  B<-princomp$scores
  
  mod_km<-NbClust(data = B,  method = "kmeans",min.nc = 2, max.nc = num_cluster_seq[length(num_cluster_seq)])
  mod_hc<-NbClust(data = B,  method = "ward.D2",min.nc = 2, max.nc = num_cluster_seq[length(num_cluster_seq)])
  
  
  mod_km_kopt<-(mod_km$Best.nc[1,])[!mod_km$Best.nc[1,]==0]
  mod_hc_kopt<-(mod_hc$Best.nc[1,])[!mod_hc$Best.nc[1,]==0]
  par(mfrow=c(1,2))
  hist(mod_km_kopt,breaks = 0:num_cluster_seq[length(num_cluster_seq)])
  hist(mod_hc_kopt,breaks = 0:num_cluster_seq[length(num_cluster_seq)])
  K_opt_km<-as.numeric(names(sort(table(mod_km_kopt),decreasing=TRUE))[1])
  K_opt_hc<-as.numeric(names(sort(table(mod_hc_kopt),decreasing=TRUE))[1])
  
  
  ciao <- mclustBIC(B,G=c(1,num_cluster_seq))
  plot(ciao)
  
  model<-Mclust(B,G=c(1,num_cluster_seq))
  
  mod_opt_hc<-hcut(B,k =K_opt_hc )
  mod_opt_km<-kmeans(B,centers =K_opt_km )
  
  mod_opt_mod<-model
  
  mean_fd<-list()
  
  mean_fd_hc<-fd(sapply(1:K_opt_hc, function(ll)colMeans(B[mod_opt_hc$cluster==ll,])%*%t(princomp$harmonics$coefs)),princomp$harmonics$basis) + 
    fd(matrix(mean_fd_ini$coefs, nrow = nrow(mean_fd_ini$coefs), ncol = K_opt_hc), mean_fd_ini$basis)
  mean_fd_km<-fd(t(mod_opt_km$centers%*%t(princomp$harmonics$coefs)),princomp$harmonics$basis) + 
    fd(matrix(mean_fd_ini$coefs, nrow = nrow(mean_fd_ini$coefs), ncol = K_opt_km), mean_fd_ini$basis)
  mean_fd_model<-fd(t(t(mod_opt_mod$parameters$mean)%*%t(princomp$harmonics$coefs)),princomp$harmonics$basis) + 
    fd(matrix(mean_fd_ini$coefs, nrow = nrow(mean_fd_ini$coefs), ncol = model$G), mean_fd_ini$basis)
  
  
  out<-list(mod_opt=list(ind_hc = mod_opt_hc,ind_km = mod_opt_km,mod_opt = mod_opt_mod),
            mean_fd=list(hc = mean_fd_hc,km = mean_fd_km,model = mean_fd_model),
            mod_km=mod_km,
            mod_hc=mod_hc)
  
}

fil_fpca_ss<-function(data,num_clusters=3:4,per_comp=0.9,grid=NULL,met="sil",...){
  
  
  basis<- create.bspline.basis(c(0,1), nbasis=20)
  loglam         = seq(-10, -4, 0.25)
  Gcvsave        = numeric()
  for(i in 1:length(loglam)){
    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    Sm.i    = smooth.basis(data$grid, data$X, fdPari)
    Gcvsave[i] = sum(Sm.i$gcv)
    
  }
  
  
  lambda_s=10^loglam[which.min(Gcvsave)]
  
  plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
       ylab=expression(GCV(lambda)), lwd=2 )
  
  fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
  X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
  
  
  sd_fd<-sd.fd(X_fd)
  mean_fd<-mean.fd(X_fd)
  X_std<-center.fd(X_fd)*(sd_fd)^-1
  princomp<-pca.fd(X_std,min(X_std$basis$nbasis, 30))
  cum<-cumsum(princomp$values)/(X_fd$basis$rangeval[2]-X_fd$basis$rangeval[1])
  num_com<-min(which(cum>per_comp))
  princomp<-pca.fd(X_std,num_com)
  B<-princomp$scores
  mod<-clValid(B, num_clusters, clMethods=c("hierarchical","kmeans"),
               validation="internal",method = "ward",
               maxitems = ncol(data$X))
  
  
  model<-Mclust(B,G=num_clusters)
  
  con<-mod@measures[1,,]
  dunn<-mod@measures[2,,]
  sil<-mod@measures[3,,]
  K_opt_con<-apply(con, 2, which.min)
  K_opt_dunn<-apply(dunn, 2, which.max)
  K_opt_sil<-apply(sil, 2, which.max)
  K_opt<-list(K_opt_con,K_opt_dunn,K_opt_sil)
  mod_opt_hc<-lapply(1:3,function(ii)cutree(mod@clusterObjs$hierarchical,num_clusters[K_opt[[ii]][1]]))
  mod_opt_km<-lapply(1:3,function(ii)mod@clusterObjs$kmeans[[K_opt[[ii]][2]]])
  mod_opt_mod<-model
  
  
  mean_fd_hc<-lapply(1:3,function(ii)fd(sapply(1:length(unique(mod_opt_hc[[ii]])), function(ll)colMeans(B[mod_opt_hc[[ii]]==ll,])%*%t(princomp$harmonics$coefs)),princomp$harmonics$basis)*sd_fd+mean_fd)
  mean_fd_km<-lapply(1:3,function(ii)fd(t(mod_opt_km[[ii]]$centers%*%t(princomp$harmonics$coefs)),princomp$harmonics$basis)*sd_fd+mean_fd)
  mean_fd_model<-fd(t(t(mod_opt_mod$parameters$mean)%*%t(princomp$harmonics$coefs)),princomp$harmonics$basis)*sd_fd+mean_fd
  
  if(met=="sil"){
    ind=3
    out<-list(mod_opt_sil=list(mod_opt_hc[[ind]],mod_opt_km[[ind]],mod_opt_mod),
              mean_fd=list(mean_fd_hc[[ind]],mean_fd_km[[ind]],mean_fd_model),
              sil=sil)
  }
  else if(met=="con"){
    ind=1
    out<-list(mod_opt_con=list(mod_opt_hc[[ind]],mod_opt_km[[ind]],mod_opt_mod),
              mean_fd=list(mean_fd_hc[[ind]],mean_fd_km[[ind]],mean_fd_model),
              con=con)
  }
  else if(met=="dunn"){
    ind=2
    out<-list(mod_opt_dunn=list(mod_opt_hc[[ind]],mod_opt_km[[ind]],mod_opt_mod),
              mean_fd=list(mean_fd_hc[[ind]],mean_fd_km[[ind]],mean_fd_model),
              dunn=dunn)
  }
  else{
    out<-list(mod_opt_con=mod_opt_con,
              mod_opt_dunn=mod_opt_dunn,
              mod_opt_sil=mod_opt_sil,
              mean_fd_hc=mean_fd_hc,
              mean_fd_km=mean_fd_km,
              mean_fd_model=mean_fd_model,
              con=con,
              dunn=dunn,
              sil=sil)
  }
  
}


fil_fpca_ms<-function(data,num_cluster_seq=2:5,per_comp_vec=0.9,met="sil",...){
  
  if(length(per_comp_vec)==1){
    mod_opt<-fil_fpca_ss(data,num_clusters = num_cluster_seq,per_comp = per_comp_vec,...)
    par_opt_vec=mod_opt[[2]]
  }
  else{
    mod<-list()
    for (ii in 1:length(per_comp_vec)) {
      mod[[ii]]<-fil_fpca_ss(data,num_clusters = num_cluster_seq,per_comp = per_comp_vec[ii],met=met)
    }
    
    par_opt_vec<-sapply(mod, "[[", 3)
    if(met=="sil"|met=="dunn")
      mod_opt<-mod[[which.max(apply(par_opt_vec,2,max))]]
    else if(met=="con")
      mod_opt<-mod[[which.min(apply(par_opt_vec,2,min))]]
    else
      mod_opt=mod
    
  }
  out<-list(mod_opt=mod_opt,par_opt_vec=par_opt_vec)
  return(out)
}

# Raw ------------------------------------------------------


raw_ms_nbclust<-function(data,num_cluster_seq=2:5){
  
  
  B<-t(data$X)
  
  mod_km<-NbClust(data = B,  method = "kmeans",min.nc = 2, max.nc = num_cluster_seq[length(num_cluster_seq)])
  mod_hc<-NbClust(data = B,  method = "ward.D2",min.nc = 2, max.nc = num_cluster_seq[length(num_cluster_seq)])
  
  
  mod_km_kopt<-(mod_km$Best.nc[1,])[!mod_km$Best.nc[1,]==0]
  mod_hc_kopt<-(mod_hc$Best.nc[1,])[!mod_hc$Best.nc[1,]==0]
  par(mfrow=c(1,2))
  hist(mod_km_kopt,breaks = 0:num_cluster_seq[length(num_cluster_seq)])
  hist(mod_hc_kopt,breaks = 0:num_cluster_seq[length(num_cluster_seq)])
  K_opt_km<-as.numeric(names(sort(table(mod_km_kopt),decreasing=TRUE))[1])
  K_opt_hc<-as.numeric(names(sort(table(mod_hc_kopt),decreasing=TRUE))[1])
  
  
  ciao <- mclustBIC(B,G=c(1,num_cluster_seq))
  plot(ciao)
  
  model<-Mclust(B,G=c(1,num_cluster_seq))
  
  mod_opt_hc<-hcut(B,k =K_opt_hc )
  mod_opt_km<-kmeans(B,centers =K_opt_km )
  
  mod_opt_mod<-model
  
  
  basis<-create.bspline.basis(c(0,1),nbasis = ncol(B),norder = 2)
  mean_fd<-list()
  mean_fd_hc<-fd(sapply(1:K_opt_hc, function(ll)colMeans(B[mod_opt_hc$cluster==ll,])),basis)
  mean_fd_km<-fd(t(mod_opt_km$centers),basis)
  mean_fd_model<-fd(mod_opt_mod$parameters$mean,basis)
  
  
  out<-list(mod_opt=list(ind_hc = mod_opt_hc,ind_km = mod_opt_km,mod_opt = mod_opt_mod),
            mean_fd=list(hc = mean_fd_hc,km = mean_fd_km,model = mean_fd_model),
            mod_km=mod_km,
            mod_hc=mod_hc)
  
  
  
  
  return(out)
  
}

raw_ms<-function(data,num_cluster_seq=2:5,met="sil"){
  
  
  B<-t(data$X)
  mod<-clValid(B, num_cluster_seq, clMethods=c("hierarchical","kmeans"),
               validation="internal",method = "ward",neighbSize = 30,
               maxitems = ncol(data$X))
  model<-Mclust(B,G=num_cluster_seq)
  
  con<-mod@measures[1,,]
  dunn<-mod@measures[2,,]
  sil<-mod@measures[3,,]
  K_opt_con<-apply(con, 2, which.min)
  K_opt_dunn<-apply(dunn, 2, which.max)
  K_opt_sil<-apply(sil, 2, which.max)
  K_opt<-list(K_opt_con,K_opt_dunn,K_opt_sil)
  mod_opt_hc<-lapply(1:3,function(ii)cutree(mod@clusterObjs$hierarchical,num_cluster_seq[K_opt[[ii]][1]]))
  mod_opt_km<-lapply(1:3,function(ii)mod@clusterObjs$kmeans[[K_opt[[ii]][2]]])
  mod_opt_mod<-model
  
  basis<-create.bspline.basis(c(0,1),nbasis = ncol(B),norder = 2)
  mean_fd<-list()
  mean_fd_hc<-lapply(1:3,function(ii)fd(sapply(1:length(unique(mod_opt_hc[[ii]])), function(ll)colMeans(B[mod_opt_hc[[ii]]==ll,])),basis))
  mean_fd_km<-lapply(1:3,function(ii)fd(t(mod_opt_km[[ii]]$centers),basis))
  mean_fd_model<-fd(mod_opt_mod$parameters$mean,basis)
  
  if(met=="sil"){
    ind=3
    out<-list(mod_opt_sil=list(mod_opt_hc[[ind]],mod_opt_km[[ind]],mod_opt_mod),
              mean_fd=list(mean_fd_hc[[ind]],mean_fd_km[[ind]],mean_fd_model),
              sil=sil)
  }
  else if(met=="con"){
    ind=1
    out<-list(mod_opt_con=list(mod_opt_hc[[ind]],mod_opt_km[[ind]],mod_opt_mod),
              mean_fd=list(mean_fd_hc[[ind]],mean_fd_km[[ind]],mean_fd_model),
              con=con)
  }
  else if(met=="dunn"){
    ind=2
    out<-list(mod_opt_dunn=list(mod_opt_hc[[ind]],mod_opt_km[[ind]],mod_opt_mod),
              mean_fd=list(mean_fd_hc[[ind]],mean_fd_km[[ind]],mean_fd_model),
              dunn=dunn)
  }
  else{
    out<-list(mod_opt_con=mod_opt_con,
              mod_opt_dunn=mod_opt_dunn,
              mod_opt_sil=mod_opt_sil,
              mean_fd_hc=mean_fd_hc,
              mean_fd_km=mean_fd_km,
              mean_fd_model=mean_fd_model,
              con=con,
              dunn=dunn,
              sil=sil)
  }
  return(out)
  
}





# Varie -------------------------------------------------------------------

loglik <- function(parameters, data, vars, FullS,W=NA,AW_vec=NA,P_tot=NA,lambda_s=NA,lambda_l=NA,CK=NA,CLC=FALSE){
  
  gamma <- vars$gamma
  gcov <- vars$gcov
  curve <- data$curve
  pi <- parameters$pi
  
  S <- FullS[data$timeindex,  ]
  
  N<-length(unique(data$curve))
  K <- dim(vars$gamma)[2]
  q <- dim(vars$gamma)[3]
  
  Gamma <- parameters$Gamma
  if(is.null(parameters$mu))mu<-t(matrix(parameters$lambda.zero,q,K)+parameters$Lambda%*%t(parameters$alpha))
  else mu<-parameters$mu
  
  loglk <- 0
  
  for(i in 1:N){
    y <- data$x[data$curve == i]
    Si <- S[data$curve == i,  ]
    ni <- dim(Si)[1]
    invvar <- diag(1/rep(parameters$sigma, ni))
    covx <- Si %*% Gamma %*% t(Si) +solve(invvar)
    # sum(sapply(1:K,function(ll)pi[ll]*(2*base::pi)^(-ni/2)*det(covx)^(-1/2)*exp(-(1/2)*t(y-Si%*%t(parameters$mu)[,ll])%*%solve(covx)%*%(y-Si%*%t(parameters$mu)[,ll]))))
    temp<-sapply(1:K,function(ll)pi[ll]*mvtnorm::dmvnorm(t(y), mean=Si%*%t(mu)[,ll], covx))
    loglk=loglk+log(sum(temp))
    
  }
  
  if(!is.na(W)){
    p_l=lambda_l*t(AW_vec)%*%abs(P_tot%*%vec(t(parameters$mu)))
    p_s=lambda_s*sum(sapply(1:K,function(ll)t(parameters$mu)[,ll]%*%W%*%t(parameters$mu)[,ll]))
    p_pi=CK*sum(sapply(1:K,function(ll)log(pi[ll])))
    # print(paste(p_s," ",p_l, " ",p_pi))
    
    ploglk<-loglk-p_l-p_s+p_pi
    if(CLC==TRUE){
      EN<--sum(sapply(1:K,function(ii)sum(sapply(1:N,function(ll)vars$piigivej[ll,ii]*log(vars$piigivej[ll,ii]+10^-200)))))
      loglk=-2*loglk+2*EN
    }
    out<-round(c(loglk,ploglk[1,1]),digits = 2)
    
    
  }
  else{
    out<-loglk
  }
  return(out)
}

