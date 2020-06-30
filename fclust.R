fitfclust <-
  function(x=NULL,curve=NULL,timeindex=NULL,data=NULL, q = 5, h = 1, p = 5,
           K = 2, tol = 0.001, maxit = 20, pert =  
             0.01, grid = seq(0, 1, length = 100), hard = F, plot= F,trace=F)
  {
    # This is the main function to implement the FClust procedure.
    library(splines)
    if (is.null(data))
      data <- list(x=x,curve=curve,timeindex=timeindex)
    initfit <- fclustinit(data = data, pert = pert, grid = grid, h = h,
                          p = p, q = q, K = K)
    # Initialize the parameters
    parameters <- initfit$parameters
    vars <- initfit$vars
    S <- initfit$S
    FullS <- initfit$FullS
    sigma.old <- 0
    sigma.new <- 1
    ind <- 1
    # Main loop. Iterates between M and E steps and stops when
    # sigma  has converged.
    while(abs(sigma.old - sigma.new)/sigma.new > tol & (ind <= maxit)) {
      parameters <- fclustMstep(parameters, data, vars, S, tol, p, hard)
      vars <- fclustEstep(parameters, data, vars, S, hard)
      sigma.old <- sigma.new
      sigma.new <- parameters$sigma[1]
      if (trace)
        print(paste("Iteration", ind,": Sigma = ",sigma.new))
      # Plot cluster mean curves.  
      if(plot){
        cluster.mean <- matrix(0,K,dim(FullS)[1])
        for(k in 1:K)
          cluster.mean[k,] <- FullS %*% (parameters$lambda.zero + parameters$Lambda %*% parameters$alpha[k,  ])
        plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
        for(k in 1:K)
          lines(grid,cluster.mean[k,], col = 4, lwd = 2)
      }
      ind <- ind + 1
    }
    # Enforce parameter constraint (7)
    cfit <- fclustconst(data, parameters, vars, FullS)
    list(data=data,parameters = cfit$parameters, vars = cfit$vars, FullS = FullS,grid=grid)
  }

fclustinit <-
  function(data, pert = 0, grid = seq(0.01, 1, length = 100), h = 1, p = 2, q = 5,K = K){
    S <- FullS <- NULL
    # This function initializes all the parameters.
    # Produce spline basis matrix
    FullS <- cbind(1, ns(grid, df = (q - 1)))
    FullS <- svd(FullS)$u
    S <- FullS[data$timeindex,  ]
    N <- length(table(data$curve))
    # Compute initial estimate of basis coefficients.
    points <- matrix(0,N,sum(q))
    for (i in 1:N){
      Si <- S[data$curve==i,]
      xi <- data$x[data$curve==i]
      points[i,] <- solve(t(Si) %*% Si + pert * diag(q)) %*% t(Si) %*%xi}
    # Use k-means to get initial cluster memberships from points.
    if(K > 1)
      class <- kmeans(points, K, 10)$cluster
    else class <- rep(1, N)
    # Initial estimates for the posterior probs of cluster membership.
    piigivej <- matrix(0, N, K)
    piigivej[col(piigivej) == class] <- 1
    # Calculate coefficeints for cluster means.
    classmean <- matrix(0,K,q)
    for (k in 1:K)
      classmean[k,] <- apply(points[class==k,],2,mean)
    # Initialize lambda.zero, Lambda and alpha as defined in paper.
    lambda.zero <- apply(classmean, 2, mean)
    Lambda <- as.matrix(svd(scale(classmean, scale = F))$v[, 1:h])
    alpha <- scale(classmean, scale = F)%*%Lambda
    # Calculate estimates for gamma and gprod.
    gamma <- t(t(points) - lambda.zero - (Lambda %*% t(alpha[class,  ])))
    gprod <- NULL
    for(i in 1:N)
      gprod <- cbind(gprod, (gamma[i,  ]) %*% t(gamma[i,  ]))
    gamma <- array(gamma[, rep(1:sum(q), rep(K, sum(q)))], c(N, K, sum(q)))
    gcov <- matrix(0, sum(q), N * sum(q))
    list(S = S, FullS = FullS, parameters = list(lambda.zero = lambda.zero,
                                                 Lambda = Lambda, alpha = alpha), vars = list(gamma = gamma,
                                                                                              gprod = gprod, gcov = gcov, piigivej = piigivej))
  }

"fclustMstep" <-
  function(parameters, data, vars, S, tol, p, hard)
  {
    # This function implements the M step of the EM algorithm.
    K <- dim(parameters$alpha)[1]
    alpha <- parameters$alpha
    Lambda <- parameters$Lambda
    gamma <- vars$gamma
    gcov <- vars$gcov
    curve <- data$curve
    piigivej <- vars$piigivej
    N <- dim(gamma)[1]
    K <- dim(alpha)[1]
    h <- dim(alpha)[2]
    n <- length(curve)
    q <- dim(S)[2]
    sigma.old <- 2
    sigma <- 1
    # Compute pi.
    if(hard)
      parameters$pi <- rep(1/K, K)
    else parameters$pi <- apply(vars$piigivej, 2, mean)
    # Compute rank p estimate of Gamma
    ind <- matrix(rep(c(rep(c(1, rep(0, q - 1)), N), 0), q)[1:(N*q^2)], N * q, q)
    gsvd <- svd(vars$gprod %*% ind/N)
    gsvd$d[ - (1:p)] <- 0
    parameters$Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
    # This loop iteratively calculates alpha and then Lambda and stops
    # when they have converged.
    loopnumber <- 1
    while((abs(sigma.old[1] - sigma[1])/sigma[1] > tol) & (loopnumber <10)){
      sigma.old <- sigma
      # Calculate lambda.zero.
      gamma.pi <- diag(S %*% t(apply(gamma * as.vector(piigivej),
                                     c(1, 3), sum)[curve,  ]))
      alpha.pi <- t(matrix(apply(alpha, 2, function(x, piigivej,K)
      {rep(1, K) %*% (piigivej * x)}
      , t(piigivej), K), N, h)[curve,  ])
      lambda.zero <- solve(t(S) %*% S) %*% t(S) %*% (data$x - diag(
        S %*% Lambda %*% alpha.pi) - gamma.pi)
      x <- data$x - S %*% lambda.zero
      # Calculate alpha.
      for(i in 1.:K) {
        S.Lam <- S %*% Lambda
        S.Lam.pi <- S.Lam * piigivej[curve, i]
        if(sum(piigivej[, i]) > 10^(-4))
          alpha[i,  ] <- solve(t(S.Lam.pi) %*% S.Lam) %*%
          t(S.Lam.pi) %*% (x - diag(S %*% t(gamma[curve, i,  ])))
        else print("Warning: empty cluster")
      }
      # Calculate Lambda given alpha. This is done by iterating
      # through each column of Lambda holding the other columns fixed.
      for(m in 1:h) {
        pi.alphasq <- apply(t(piigivej) * (alpha^2)[, m], 2,sum)[curve]
        pi.alpha <- apply(t(piigivej) * alpha[, m], 2, sum)[curve]
        S.Lambda <- t(S %*% Lambda)
        if(h != 1) {
          temp <- NULL
          for(i in 1:K) {
            temp <- cbind(temp, as.vector(rep(1, h - 1) %*% matrix((S.Lambda *
                                                                      alpha[i,  ])[ - m,  ], h - 1,dim(S)[1])) * alpha[i, m])
          }
          otherLambda <- (temp * piigivej[curve,  ])%*%rep(1, K)
        }
        else otherLambda <- 0
        gamma.pi.alpha <- apply(gamma * as.vector(piigivej) *
                                  rep(alpha[, m], rep(N, K)), c(1, 3), sum)[curve,  ]
        Lambda[, m] <- solve(t(S * pi.alphasq) %*% S) %*% t(S) %*%
          (x * pi.alpha - otherLambda - (S *gamma.pi.alpha) %*% rep(1, sum(q)))
      }
      # Calculate sigma 
      ind <- matrix(rep(c(rep(F, q), rep(T, N * q)),N)
                    [1:(N * N * q)], N, N * q, byrow = T)[rep(1:N, table(curve)),]
      mat1 <- matrix(rep(S, N), n, N * q)
      mat3 <- t(mat1)
      mat3[t(ind)] <- 0
      ind2 <- matrix(rep(c(rep(F, q), rep(T, N * q)),
                         N)[1:(N * N * q)], N, N * q, byrow = T)[rep(1:N, rep(q, N)),]
      mat2 <- matrix(rep(t(gcov), N), N * q, N * q,byrow = T)
      mat2[ind2] <- 0
      econtrib <- 0
      for(i2 in 1:K) {
        vect1 <- x - S %*% Lambda %*% alpha[
          i2,  ] - (S * gamma[curve, i2,  ]) %*% rep(1, q)
        econtrib <- econtrib + t(piigivej[curve,i2] * vect1) %*% vect1
      }
      sigma <- as.vector((econtrib + sum(diag(mat1 %*% mat2 %*% mat3)))/n)
      loopnumber <- loopnumber + 1
    }
    parameters$lambda.zero <- as.vector(lambda.zero)
    parameters$alpha <- alpha
    parameters$Lambda <- Lambda
    parameters$sigma <- sigma
    parameters
  }

"fclustEstep" <-
  function(parameters, data, vars, S, hard)
  {
    # This function performs the E step of the EM algorithm by
    # calculating the expected values of gamma and gamma %*% t(gamma)
    # given the current parameter estimates.
    par <- parameters
    N <- dim(vars$gamma)[1]
    K <- dim(vars$gamma)[2]
    q <- dim(vars$gamma)[3]
    Gamma <- par$Gamma
    Lambda.alpha <- par$lambda.zero + par$Lambda %*% t(par$alpha)
    vars$gprod <- vars$gcov <- NULL
    for(j in 1:N) {
      # Calculate expected value of gamma.
      Sj <- S[data$curve == j,  ]
      nj <- sum(data$curve == j)
      invvar <- diag(1/rep(par$sigma, nj))
      Cgamma <- Gamma - Gamma %*% t(Sj) %*% solve(diag(nj) + Sj %*%
                                                    Gamma %*% t(Sj) %*% invvar) %*% invvar %*% Sj %*% Gamma
      centx <- data$x[data$curve == j] - Sj %*% Lambda.alpha
      vars$gamma[j,  ,  ] <- t(Cgamma %*% t(Sj) %*% invvar %*% centx)
      # Calculate pi i given j.
      covx <- Sj %*% par$Gamma %*% t(Sj) + solve(invvar)
      d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * par$pi
      vars$piigivej[j,  ] <- d/sum(d)
      if(hard) {
        m <- order( - d)[1]
        vars$piigivej[j,  ] <- 0
        vars$piigivej[j, m] <- 1
      }
      # Calculate expected value of gamma %*% t(gamma).
      vars$gprod <- cbind(vars$gprod, t(matrix(vars$gamma[j,  ,  ],
                                               K, q)) %*% (matrix(vars$gamma[j,  ,  ], K, q) * vars$
                                                             piigivej[j,  ]) + Cgamma)
      vars$gcov <- cbind(vars$gcov, Cgamma)
    }
    vars
  }


"fclustconst" <-
  function(data, parameters, vars, S)
  {
    # This function enforces the constraint (7) from the paper on the
    # parameters. This means that the alphas can be interpreted as the
    # number of standard deviations that the groups are apart etc.
    par <- parameters
    A <- t(S) %*% solve(par$sigma * diag(dim(S)[1]) + S %*% par$Gamma %*%
                          t(S)) %*% S
    svdA <- svd(A)
    sqrtA <- diag(sqrt(svdA$d)) %*% t(svdA$u)
    negsqrtA <- svdA$u %*% diag(1/sqrt(svdA$d))
    finalsvd <- svd(sqrtA %*% par$Lambda)
    par$Lambda <- negsqrtA %*% finalsvd$u
    if(dim(par$Lambda)[2] > 1)
      par$alpha <- t(diag(finalsvd$d) %*% t(finalsvd$v) %*% t(par$alpha))
    else par$alpha <- t(finalsvd$d * t(finalsvd$v) %*% t(par$alpha))
    meanalpha <- apply(par$alpha, 2, mean)
    par$alpha <- t(t(par$alpha) - meanalpha)
    par$lambda.zero <- par$lambda.zero + par$Lambda %*% meanalpha
    list(parameters = par, vars = vars)
  }

"nummax" <-
  function(X)
  {
    ind <- rep(1, dim(X)[1])
    m <- X[, 1]
    if(dim(X)[2] > 1)
      for(i in 2:dim(X)[2]) {
        test <- X[, i] > m
        ind[test] <- i
        m[test] <- X[test, i]
      }
    list(ind = ind, max = m)
  }

"fclust.pred" <-
  function(fit,data=NULL,reweight=F)
  {
    # This function produces the alpha hats used to provide a low
    # dimensional pictorial respresentation of each curve. It also
    # produces a class prediction for each curve. It takes as
    # input the fit from fldafit (for predictions on the original data)
    # or the fit and a new data set (for predictions on new data).
    if (is.null(data))
      data <- fit$data
    FullS <- fit$FullS
    par <- fit$parameters
    curve <- data$curve
    time <- data$time
    N <- length(table(curve))
    h <- dim(par$alpha)[2]
    alpha.hat <- matrix(0, N, h)
    K <- dim(fit$par$alpha)[1]
    distance <- matrix(0, N, K)
    Calpha <- array(0, c(N, h, h))
    for(i in 1:N) {
      Sij <- FullS[time[curve == i],  ]
      xij <- data$x[curve == i]
      n <- length(xij)
      Sigma <- par$sigma * diag(n) + Sij %*% par$Gamma %*% t(Sij)
      # Calculate covariance for each alpha hat.
      InvCalpha <- t(par$Lambda) %*% t(Sij) %*% solve(Sigma) %*% Sij %*%
        par$Lambda
      Calpha[i,  ,  ] <- solve(InvCalpha)
      # Calculate each alpha hat.
      alpha.hat[i,  ] <- Calpha[i,  ,  ] %*% t(par$Lambda) %*% t(
        Sij) %*% solve(Sigma) %*% (xij - Sij %*% par$lambda.zero)
      # Calculate the matrix of distances, relative to the
      # appropriate metric of each curve from each class centroid. 
      for (k in 1:K){
        y <- as.vector(alpha.hat[i,])-fit$par$alpha[k,]
        distance[i,k] <- t(y)%*%InvCalpha %*%y}}
    # Calculate final class predictions for each curve.
    class.pred <- rep(1, N)
    log.pi <- log(fit$par$pi)
    if (!reweight)
      log.pi <- rep(0,K)
    probs <- t(exp(log.pi-t(distance)/2))
    probs <- probs/apply(probs,1,sum)
    m <- probs[,1]
    if(K != 1)
      for(k in 2:K) {
        test <- (probs[, k] > m)
        class.pred[test] <- k
        m[test] <- probs[test, k]
      }
    list(Calpha = Calpha, alpha.hat = alpha.hat, class.pred = class.pred,
         distance = distance, m = m,probs=probs)
  }

"fclust.curvepred" <-
  function(fit, data=NULL, index=NULL, tau = 0.95, tau1 = 0.975)
  {
    if (is.null(data))
      data <- fit$data
    if (is.null(index))
      index <- 1:length(table(data$curve))
    tau2 <- tau/tau1
    sigma <- fit$par$sigma
    Gamma <- fit$par$Gamma
    Lambda <- fit$par$Lambda
    alpha <- fit$par$alpha
    lambda.zero <- as.vector(fit$par$lambda.zero)
    S <- fit$FullS
    N <- length(index)
    upci <-lowci <- uppi <- lowpi <- gpred <- matrix(0,N,nrow(S))
    etapred <- matrix(0,N,ncol(S))
    ind <- 1
    Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
    for (i in index){
      y <- data$x[data$curve == i]
      Si <- S[data$time[data$curve == i],  ]
      ni <- dim(Si)[1]
      invvar <- diag(1/rep(sigma, ni))
      covx <- Si %*% Gamma %*% t(Si) + solve(invvar)
      centx <- data$x[data$curve == i] - Si %*% Lambda.alpha
      d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * fit$par$pi
      pi <- d/sum(d)
      K <- length(pi)
      mu <- lambda.zero + Lambda %*% t(alpha * pi) %*% rep(1, K)
      cov <- (Gamma - Gamma %*% t(Si) %*% solve(diag(ni) + Si %*% Gamma %*%
                                                  t(Si)/sigma) %*% Si %*% Gamma/sigma)/sigma
      etapred[ind,] <- mu + cov %*% t(Si) %*% (y - Si %*% mu)
      ord <- order( - pi)
      numb <- sum(cumsum(pi[ord]) <= tau1) + 1
      v <- diag(S %*% (cov * sigma) %*% t(S))
      pse <- sqrt(v + sigma)
      se <- sqrt(v)
      lower.p <- upper.p <- lower.c <- upper.c <- matrix(0, nrow(S), numb)
      for(j in 1:numb) {
        mu <- lambda.zero + Lambda %*% alpha[ord[j],  ]
        mean <- S %*% (mu + cov %*% t(Si) %*% (y - Si %*% mu))
        upper.p[, j] <- mean + qnorm(tau2) * pse
        lower.p[, j] <- mean - qnorm(tau2) * pse
        upper.c[, j] <- mean + qnorm(tau2) * se
        lower.c[, j] <- mean - qnorm(tau2) * se
      }
      upci[ind,] <- nummax(upper.c)$max
      lowci[ind,] <-  - nummax( - lower.c)$max
      uppi[ind,] <- nummax(upper.p)$max
      lowpi[ind,] <-  - nummax( - lower.p)$max
      gpred[ind,] <- as.vector(S %*%etapred[ind,])
      ind <- ind+1
    }
    meancurves <- S%*%Lambda.alpha
    list(etapred = etapred, gpred = gpred,  upci = upci,lowci = lowci,  uppi = uppi, lowpi = lowpi,index=index,grid=fit$grid,data=data,meancurves=meancurves)
  }

"fclust.discrim" <-
  function(fit,absvalue=F){
    S <- fit$FullS
    sigma <- fit$par$sigma
    n <- nrow(S)
    Gamma <- fit$par$Gamma
    Sigma <- S%*%Gamma%*%t(S)+sigma*diag(n)
    Lambda <- fit$par$Lambda
    discrim <- solve(Sigma)%*%S%*%Lambda
    if (absvalue)
      discrim <- abs(discrim)
    n <- ncol(discrim)
    nrows <- ceiling(sqrt(n))
    par(mfrow=c(nrows,nrows))
    for (i in 1:n){
      plot(fit$grid,discrim[,i],ylab=paste("Discriminant Function ",i),xlab="Time",type='n')
      lines(fit$grid,discrim[,i],lwd=3)
      abline(0,0)}}

"fclust.plotcurves" <-
  function(object=NULL,fit=NULL,index=NULL,ci=T,pi=T,clustermean=F){
    if (is.null(object))
      object <- fclust.curvepred(fit)
    if (is.null(index))
      index <- 1:length(table(object$data$curve))
    r <- ceiling(sqrt(length(index)))
    par(mfrow=c(r,r))
    for (i in index){
      grid <- object$grid
      upci <- object$upci[i,]
      uppi <- object$uppi[i,]
      lowci <- object$lowci[i,]
      lowpi <- object$lowpi[i,]
      gpred <- object$gpred[i,]
      meancurves <- (object$mean)
      if (clustermean)
        yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
      else
        yrange <- c(min(lowpi),max(uppi))
      plot(grid,grid,ylim=yrange,ylab="Predictions",xlab="Time",type='n',
           main=paste("Curve ",i))
      if (clustermean)
        for (k  in 1:ncol(meancurves))
          lines(grid,meancurves[,k],col=6,lty=2,lwd=2)
      if (ci){
        lines(grid,upci,col=3)
        lines(grid,lowci,col=3)}
      if (pi){
        lines(grid,uppi,col=4)
        lines(grid,lowpi,col=4)}
      lines(grid,gpred,col=2,lwd=2)
      lines(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],lwd=2)
      points(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],pch=19,cex=1.5)
    }
  }

"simdata" <-
  structure(list(x = c(-0.227117840025744, -0.195701325120710, 
                       -0.123029609099774, -0.118619923345801, -0.110854003404887, -0.106908797785542, 
                       -0.107349467323043, -0.136947437978017, 0.147316597820139, 0.342213089281111, 
                       -0.374637137339911, -0.132569470675019, -0.106261313334034, -0.184045390440344, 
                       -0.176687325709244, -0.205370648502014, -0.207503074221494, -0.223960531903988, 
                       -0.179282627425149, 0.228839712560303, -0.406436079640189, -0.252485300162310, 
                       -0.130893822115844, -0.119483603887501, -0.157911812412119, -0.191173118294793, 
                       -0.186391707679505, -0.166101182172451, 0.017989520648983, 0.102677621538723, 
                       -0.329746791407803, -0.280156068523990, -0.194771696246444, -0.149493252880510, 
                       -0.171596928524117, -0.197948023475606, -0.217010552084032, -0.131802072904116, 
                       0.256737323621704, 0.358579862077912, -0.226652568726055, -0.194959419590792, 
                       -0.143850156023500, -0.100379181559697, -0.106534204920563, -0.182279005060978, 
                       -0.191699200147385, -0.209233185188917, -0.131366576347318, 0.168550341932226, 
                       -0.313790919034796, -0.116777642865946, -0.101128217858913, -0.104795564464629, 
                       -0.094204404109515, -0.106849583723510, -0.114156343156498, -0.113507166276034, 
                       -0.137597617774813, 0.340065783702146, -0.129486325727864, -0.130020108924852, 
                       -0.136359618769366, -0.162790169199692, -0.174869930544662, -0.193164071088511, 
                       -0.203685565671835, -0.093588403559253, 0.232631809906245, 0.465850115143711, 
                       -0.294711458343317, -0.0948725387147703, -0.110915079204378, 
                       -0.114290893388158, -0.160756593743237, -0.179658994991600, -0.191269772791207, 
                       -0.223781346510862, 0.203603299174188, 0.331811990712914, -0.467522085682574, 
                       -0.245553797585632, -0.158979981457772, -0.209709876424111, -0.220642628521095, 
                       -0.169252620230874, -0.0887299961627144, -0.0386754613677018, 
                       0.147942919391726, 0.414151222903763, -0.112767077728689, -0.124058844743678, 
                       -0.107685638343733, -0.128089652858022, -0.206697524949735, -0.152130234525033, 
                       -0.0748596280764934, -0.0450328357374213, 0.237789276941251, 
                       0.317903890876979, -0.325736649236121, -0.305421771375266, -0.136924941126976, 
                       -0.137501625477845, -0.148820202759715, -0.155879507713238, -0.196535108694717, 
                       -0.213591244573685, 0.185888155688001, 0.436176736410145, -0.238859561543531, 
                       -0.175022237353340, -0.150749556903914, -0.120811412581943, -0.0948573655562948, 
                       -0.114608637387264, -0.191172926480358, -0.179682834570726, 0.0754949932557839, 
                       0.119940540722110, -0.366315219441288, -0.178066636999900, -0.155689951246938, 
                       -0.107700926401963, -0.191011693833892, -0.229049245707688, -0.219123020434303, 
                       -0.129898949730155, 0.140321315344435, 0.478295571715595, -0.281431065900907, 
                       -0.0706397378256686, -0.106080777103587, -0.115991974742118, 
                       -0.159791068256919, -0.192988470004975, -0.139795187549111, 0.0973676688021566, 
                       0.38102556559798, 0.496720144241713, -0.211171865959491, -0.138293017933202, 
                       -0.166794027612996, -0.181076913337912, -0.182749238304548, -0.168885929237816, 
                       -0.0587760637453383, 0.0316074224492748, 0.155584908558565, 0.555522489857623, 
                       -0.290406447096274, -0.220314893140696, -0.112347777803585, -0.115267695281464, 
                       -0.160607886057124, -0.207677052224203, -0.21910912393415, -0.20626539566681, 
                       -0.0637741762261909, 0.108964406226199, -0.290948473145957, -0.155906091007454, 
                       -0.134894200821459, -0.158232356414957, -0.0414731840912619, 
                       0.035033756130331, 0.0542306741882772, 0.212790645279318, 0.397827685000681, 
                       0.51291730430605, -0.156524715515184, -0.130045132779869, -0.155544030068691, 
                       -0.155970783481934, -0.169382977727972, -0.196026213763414, -0.207036246402719, 
                       -0.205984182463623, -0.176546877232985, 0.387796222378189, -0.249795263531387, 
                       -0.185798771958970, -0.146688823530977, -0.188600638812819, -0.202593538140893, 
                       -0.210767955364829, -0.234011929530135, -0.13170472959584, 0.125620597102558, 
                       0.279205017399625, -0.399851233230859, -0.219026331078527, -0.170436950199829, 
                       -0.102921774569798, -0.127791067850398, -0.200706183728387, -0.173932253628848, 
                       -0.163008877591484, -0.120593207428063, 0.495016377680583, -0.286600479202477, 
                       -0.119262974845128, -0.124745359718768, -0.116416916732292, -0.203009989022764, 
                       -0.194560224138502, -0.171206811411954, -0.0231962559831943, 
                       0.088131052166345, 0.324252569958893, -0.0931900422255536, -0.178908326844773, 
                       -0.192336968551136, -0.206658811852580, -0.131476853218849, -0.102582260078396, 
                       -0.084498687186169, -0.0345087225881057, 0.109762905726129, 0.122130319108726, 
                       -0.454583210403097, -0.304708227418333, -0.257661666700366, -0.167608549363106, 
                       -0.120998634233420, -0.141831428444438, -0.197273334232284, -0.00197764264859539, 
                       0.0971596256025763, 0.444728583038965, -0.359705490898084, -0.256955148013682, 
                       -0.146425421489230, -0.218327769987794, -0.224592291014539, -0.238186585845544, 
                       -0.196196645095266, -0.00570295167456461, 0.0558610176444916, 
                       0.128661145304187, -0.418582569739349, -0.399575029953990, -0.328082188483067, 
                       -0.289594344484002, -0.213606099580371, -0.079966297024145, -0.156335800652327, 
                       -0.226997343996872, 0.0808948589219096, 0.531521213665854, -0.199251238106442, 
                       -0.186073436253551, -0.168433847004301, -0.128458227364381, -0.119513370328817, 
                       -0.15659155256917, -0.177241291157763, -0.192573532673733, -0.176340869411284, 
                       0.197742997138806, -0.233437597630301, -0.0768722224347693, -0.106420303349253, 
                       -0.106760502450435, -0.109720336848384, -0.127683260243216, -0.159233350259349, 
                       -0.154148969576236, -0.135234235822275, -0.0783323295016301, 
                       -0.352976275660130, -0.165281489761331, -0.106852974752975, -0.120282696258084, 
                       -0.0900729622949594, -0.108095122962697, -0.115086268548261, 
                       -0.188557742897746, -0.194469463330198, -0.162927560227355, -0.251412939971511, 
                       -0.152001981565379, -0.142398824435427, -0.118530291927102, -0.101114832748401, 
                       -0.179537684037522, -0.0825895821729416, 0.0449213623130559, 
                       0.245654552891542, 0.312217600980186, -0.378711193336546, -0.101049793095541, 
                       -0.113953924332348, -0.156151729713191, -0.223358035624415, -0.164319853879623, 
                       0.000856493745469496, 0.282248335037168, 0.308281698774851, 0.487909769779855, 
                       -0.323956437563084, -0.156880466203510, -0.113600022700249, -0.135290936426864, 
                       -0.191702652907591, -0.181490087463943, -0.165096296426790, -0.116664880059740, 
                       0.228968440123362, 0.285225778499167, -0.358106129154346, -0.247082740509186, 
                       -0.207534453544336, -0.223753386709386, -0.146350409916299, -0.176900270737067, 
                       -0.190458423537747, 0.0260398632656352, 0.134733082201608, 0.228337230203115, 
                       -0.190187320713332, -0.121251374569518, -0.117727365105491, -0.121993285934650, 
                       -0.11062489628138, -0.201116039578890, -0.0313886533690395, -0.00682049560319782, 
                       0.141526122484713, 0.330995264333811, -0.381377741444594, -0.291244001380736, 
                       -0.192009372058398, -0.162136668333722, -0.144486459566970, -0.107144638860840, 
                       -0.205388296163639, -0.0626299081559675, 0.411127357348002, 0.424261078644716, 
                       -0.261107710599625, -0.127263702045977, -0.19619577114638, -0.213093493556485, 
                       -0.194774726896583, -0.14388328046612, -0.138167916177372, 0.174560503379344, 
                       0.238645373821321, 0.464357626801585, -0.233199980491786, -0.185141893381162, 
                       -0.0857527163360715, -0.0854086316244292, -0.092761990823785, 
                       -0.0931739718623003, -0.210296065991791, -0.198063526130457, 
                       -0.0808428213414711, 0.256282450144409, -0.309022339788325, -0.166168687431289, 
                       -0.108276868230148, -0.157561151664909, -0.220907667862605, 0.116098193751447, 
                       0.271612020353032, 0.307372178178702, 0.43322553432155, 0.464702850668929, 
                       -0.35380336024507, -0.232762467143366, -0.233398498276365, -0.166782272455403, 
                       -0.157385201509253, -0.159022853485096, -0.185931911193076, -0.194953413185561, 
                       -0.0302504228154871, 0.204485358729827, -0.258401133433710, -0.145090790375331, 
                       -0.115575906598350, -0.135055923034147, -0.199044996713163, -0.176151135900161, 
                       -0.187643789271641, -0.147133289680793, 0.0326542393449244, 0.430477528181193, 
                       -0.331200291000149, -0.128722874170918, -0.109245755782248, -0.104183566734898, 
                       -0.125399534657987, -0.122202266569091, -0.184381268388992, -0.178479750089626, 
                       -0.17198079672465, -0.0930406282630965, -0.203041548955876, -0.0791067897757055, 
                       -0.0756352849424581, -0.187285535565779, -0.219649455554267, 
                       -0.204365024714019, -0.196861722813166, -0.0446299111261032, 
                       0.00406490826235925, 0.186641334906245, -0.304497226817912, -0.276608990743012, 
                       -0.252416142391013, -0.166108469942520, -0.117000338106418, -0.157451822179411, 
                       -0.154932689174077, -0.186796199638869, 0.0629520004317411, 0.327079634477786, 
                       -0.209111477685358, -0.127355775366502, -0.138003886806337, -0.106588004853118, 
                       -0.141432959810707, -0.138191763755819, -0.194943023466968, -0.208677822541162, 
                       0.0146127341972763, 0.254084214986552, -0.154937362229819, -0.107557230984147, 
                       -0.128911292441997, -0.108265595388336, -0.0972065305960794, 
                       -0.198848339248592, -0.144493203835175, -0.127644182171595, 0.0134502108900903, 
                       0.124670074165381, -0.245855478501020, -0.0998303247691579, -0.162129013306720, 
                       -0.195290395791238, -0.213403515093102, -0.187742087568606, -0.0594370281039981, 
                       0.13861977388641, 0.249122754855277, 0.35333681519369, -0.383534051446049, 
                       -0.337105891577934, -0.268760445129686, -0.223563647432947, -0.187341764032434, 
                       -0.133750457373021, -0.0990960990868987, -0.105793964717356, 
                       -0.197739313397286, -0.201583639800325, -0.403453357083668, -0.334359627979758, 
                       -0.123047612697346, -0.168076738661661, -0.196851530695217, -0.195770584738134, 
                       -0.127781545323436, -0.00138040098797243, 0.254865585978388, 
                       0.357849546372697, -0.410062967674191, -0.194539179448576, -0.162021670718332, 
                       -0.113893440894988, -0.182019633345952, -0.19681832151645, -0.167796353259233, 
                       -0.0986264620919684, 0.0145880192408087, 0.000852395835108693, 
                       -0.337756152762277, -0.288876923550474, -0.273713731979953, -0.127827303667975, 
                       -0.114208203637594, -0.139389680399877, -0.200024901975884, 0.154044561709217, 
                       0.526951764744455, 0.553860216098293, -0.467616845874366, -0.340601117677694, 
                       -0.205281183920167, -0.141800695467282, -0.226455564542709, -0.201159315944254, 
                       -0.1143799109681, -0.00098168967092759, 0.158534603947694, 0.248783978558958, 
                       0.0553287059125787, 0.0215246054422910, 0.0239613786013785, -0.0176461865447992, 
                       -0.0191688793145776, -0.00730960984810723, 0.00878934904284103, 
                       -0.0350527629503127, -0.0346358146176491, -0.0296955272586848, 
                       0.0110990699552697, -0.000952413819932563, 0.00142499831522622, 
                       -0.00342323171095836, 0.00438588540292069, -0.0223721905235238, 
                       -0.0131215092016118, 0.0129529001556294, -0.00180874179435203, 
                       0.042312234025874, 0.0143706652503577, 0.00532707892538855, -0.0049342694664682, 
                       -0.0144482297618435, -0.0240816432962334, -0.000135491435180247, 
                       -0.0113955327736495, -0.0109439317195270, 0.00951796187222357, 
                       0.0163202637078802, -0.0198606398363495, -0.020635120676532, 
                       -0.0153767479007658, -0.00486360430473615, -0.0102702756445511, 
                       -0.0289056176027155, -0.0239719675169262, -0.0172635358534988, 
                       -0.0166558349512558, -0.0298517909475893, 0.0184147010293975, 
                       -0.000463354037578969, -0.00843094145752214, 0.0259775422193117, 
                       0.028218108776461, 0.043027727565836, 0.0177876626833985, -0.00209464809102063, 
                       0.00119053268513568, 0.0369744467124425, -0.00720553911161553, 
                       -0.0151499711639636, -0.0284534208666103, -0.009622788707331, 
                       0.0161730783373770, 0.0323081379127002, 0.0259789846733687, 0.0328561802850742, 
                       0.0102371269567881, 0.0182876433615863, -0.0256289639842253, 
                       -0.0352357430051817, 0.0128692352011459, 0.0251996890570174, 
                       0.0110481162699525, 0.0107128616176025, -0.0086034745208899, 
                       -0.0357225631776305, -0.0608928376837651, 0.0629855659989359, 
                       0.0206071124605609, -0.00507257075420088, 0.00105083578850915, 
                       0.0101881936420306, -0.00294732499189758, 0.00784036852778495, 
                       0.0148645714263133, -0.00168913769964944, 0.0464950451773661, 
                       0.0371357804104615, 0.00440951807571682, 0.0063257574756019, 
                       -0.00690308496033264, -0.0127101315235538, -0.00642222391819524, 
                       -0.000953975000636369, -0.0103237791754173, -0.000157077623329786, 
                       0.0177794283768949, -0.0390410564056978, -0.036583451988664, 
                       -0.0565764259556362, -0.03357769274774, -0.0404416771728977, 
                       -0.00119215090711323, -0.00242262677533991, 0.0024386484654713, 
                       0.000518421577880747, -0.0157624852731680, 0.00297059583941147, 
                       -0.0292535990710019, -0.0261833779750313, -0.0124220376602145, 
                       0.0186116165045486, 0.0270116739543789, 0.00278702139867078, 
                       0.00461217740446334, -0.0163692493133565, -0.0134683246870467, 
                       -0.0118760549656587, -0.0143428534689524, -0.00154621883673825, 
                       0.0122181649173885, -0.00675496510840951, -0.00831456807073226, 
                       -0.0118703332287872, -0.0280348049127046, -0.0353312859539545, 
                       -0.0324108550270698, -0.0318826379901208, -0.0102942532966120, 
                       0.00716491666055001, 0.0066479365649925, 0.0263461025271677, 
                       0.0254512096720809, -0.00830943007307247, 0.00515257993210898, 
                       -0.0399028450605823, -0.0314372612152567, 0.0134045044082605, 
                       -0.0238106633743421, -0.0150826480512936, 0.0232642844997328, 
                       0.0518754025752096, 0.0589901087133986, -0.0141618913856349, 
                       -0.0180169458288145, -0.014090111280304, 0.00622379166698127, 
                       0.0072671030183475, -0.0295185439427869, -5.45286785335244e-05, 
                       -0.0270962452836659, -0.0255441235826307, -0.0235430369204862, 
                       -0.0185205480703060, -0.020367305484681, 0.0443908083782479, 
                       0.057892368385636, 0.0616258024119745, -0.00938962662262531, 
                       0.0176101729547615, -0.0245162489771747, 0.000667033018355962, 
                       -0.0340306898151468, -0.0458999905084088, -0.0235290084600596, 
                       -0.0244923823880294, 0.00818426138625731, 0.0240572963182675, 
                       -0.0538742924810529, -0.00530525005164201, 0.0149225607254, 0.0231206943231519, 
                       0.0209580945380243, -0.0409476490604760, -0.0278475671284697, 
                       -0.0348738947594903, -0.0541561685544011, -0.0500253980100138, 
                       0.108852925351926, 0.0791851403270575, 0.0164913311613404, -0.0357430088910859, 
                       0.0265907572387279, 0.0321196389260411, 0.00812212169141438, 
                       0.0231386122689169, -0.0100993273908665, 0.0148337288584504, 
                       0.0237519498904415, 0.0232655049596983, -0.00315611525870672, 
                       -0.0182097458463929, -0.0250820510435556, -0.0237785871603054, 
                       -0.0143865605785534, 0.0176929867493753, 0.00559086152356567, 
                       0.0365900938948946, 0.0725413289079942, 0.00847775220326222, 
                       0.0116704657376147, -0.00535217196916652, -0.0261805266809155, 
                       0.00966422338718637, -0.00544047151885465, -0.0314687869841454, 
                       -0.0215090737485400, 0.0023541340214244, -0.0317045838502471, 
                       0.0169934128196194, -0.00357161700195238, 8.29888951480312e-05, 
                       -0.00689143660466761, -0.0114808816859191, -0.00752078861851322, 
                       0.00881308173460774, 0.00456660066177290, 0.0202523775268966, 
                       -0.0643903701458935, -0.00132494042775946, -0.00958515505986945, 
                       0.00845340622478185, -0.00271761403286989, -0.0308326908600265, 
                       -0.0191672573737609, -0.0101604706438079, 0.0239723135820425, 
                       0.0272102422943043, 0.0056382746842277, -0.00542724129827099, 
                       -0.00936678502309076, -0.00398235776426829, 0.0342902102637111, 
                       0.0429914641123852, 0.0291034096437177, 0.0324847484276853, 0.030368576627736, 
                       0.0111513720556588, 0.00517083360824700, -0.0247401208560401, 
                       -0.0318358675335491, -0.0295561909256352, -0.0408570585805763, 
                       -0.0528791326550513, -0.0196349166496362, -0.0349120535701848, 
                       -0.0359461579261282, -0.0437609924454309, -0.00288687514902745, 
                       0.00731921993890389, 0.0167571723246924, 0.0233454350708218, 
                       0.0318169082200929, 0.0281963570688959, 0.0506079191022968, 0.0582116376591479, 
                       0.070367872540992, 0.043527023687034, -0.00649479787890526, -0.00444883239221610, 
                       0.0174604878972641, 0.00675542826279528, 0.00549051666873473, 
                       -0.0307035050700433, -0.00490307170718939, 0.00281447320544373, 
                       0.00263436860057402, 0.0126143400578613, -0.0173094942882652, 
                       -0.0212338074108356, -0.00173470796233160, 0.0147299176927094, 
                       0.0309281609372946, 0.0183195100235955, -0.00694645438425976, 
                       0.00378692310389241, 0.0152769725444509, 0.0124086734454976, 
                       0.00878378227879518, 0.0271247717484655, 0.0634609121188237, 
                       0.0490235613085486, 0.00891374041550685, 0.0168778060762436, 
                       0.0193967436522599, 0.0153342400386201, 0.0240256954404721, 0.0160073957535563, 
                       -0.0488617118112665, -0.0447079146969686, -0.0291955140364235, 
                       -0.0150909480186029, 0.0348944164725763, 0.0417604830705330, 
                       0.0374818264095243, 0.00876074962097476, 0.0139010609676386, 
                       0.00734062472377484, -0.0166397661189595, 0.0128922628121788, 
                       -0.0264481565105376, -0.0419037630823800, -0.0435989114416279, 
                       -0.0248555118960409, -0.0121380842233070, -0.0162510287969617, 
                       -0.00354036763043197, -0.0191266918530752, -0.0405654225850345, 
                       0.0279479674449364, 0.00530179776733666, -0.00434769335092900, 
                       -0.0185788147549564, -0.0235054079738696, -0.0185438578332777, 
                       -0.0178700306901307, -0.0411663572253689, -0.0434927632923373, 
                       -0.0183140126450316, -0.000563734936951124, -0.0263996143022747, 
                       -0.0305175247457801, -0.0165630736797031, -0.0266641427443859, 
                       -0.0202167859200049, -0.0153252380798637, 0.000141259634329237, 
                       0.000377706755727620, -0.0301808172611747, -0.0296942954208264, 
                       -0.0199948309462138, -0.0107166796376800, -0.030457006259209, 
                       -0.0344933486733539, -0.0213810148398748, -0.00180082254704368, 
                       -0.00587022620451189, 0.000527130826129967, 0.0435118882088142, 
                       0.000124175764797906, -0.0330382185299849, -0.00886556387605064, 
                       0.00860111128189471, -0.00353716370281205, -0.0227203563171015, 
                       -0.022602386344711, -0.0500984166847659, -0.0242015332230422, 
                       -0.0190346284541696, -0.00798568802670223, 0.00555638147758014, 
                       -0.0200720317217405, 0.00546935725826135, -0.00101071131899548, 
                       -0.0139147060420053, -0.0142389655570577, -0.0102756513886184, 
                       -0.0123288834407095, -0.0128582407814221, -0.0236995403962298, 
                       -0.0344752504955887, -0.0139979795526721, 0.00626833627078837, 
                       0.0148793542945203, 0.0368054094321584, 0.0377873237995442, 0.0196285802410993, 
                       -0.0499312376878418, -0.00216036484337487, -0.0138713673119948, 
                       -0.0174134057775548, -0.0163310184759559, 0.0058287874449919, 
                       0.00137887340193271, 0.00516545169829899, 0.00169440014622157, 
                       -0.0199444294405613, -0.0334866776837751, -0.0209787321651085, 
                       -0.0414651463159398, -0.0166678175000578, 0.00394513945530166, 
                       -0.000999504916693147, 0.0437205437903712, 0.00384724912298903, 
                       0.0133200254593556, -0.0258137689031911, -0.0352310836343037, 
                       -0.0218638464424987, 0.0245833672974279, 0.0238849367500967, 
                       0.0148482935652693, -0.0103438747229818, 0.00233285352216142, 
                       0.0123392359921639, 0.00792724769782464, 0.0146015829934379, 
                       0.00195962484270742, -0.0400900694726538, 0.0400644277243114, 
                       0.0352288116277324, 0.0267494221797731, 0.0390116347698599, -0.00598954140647742, 
                       -0.00168119191539172, -0.0164730378818075, -0.0114024889709441, 
                       -0.0122829514339268, 0.00526370750309657, 0.0196149083610763, 
                       -0.0220877620446392, 0.025333626680386, 0.000525106365000314, 
                       0.0127830667675974, 0.0418025765064845, -0.0128100293992435, 
                       -0.00498348249469088, 0.00418801643890062, -0.0250643216318743, 
                       -0.0136131242752714, -0.0158623804730166, 0.00669225391072256, 
                       0.0123125833845898, -0.00225500789742969, 0.00905849315064329, 
                       0.00157647458011261, -0.0201014528047367, -0.0408964494724469, 
                       -0.00660378745163734, 0.00242672039927597, 0.0186136872149079, 
                       0.0268361549534905, 0.00827878664532323, 0.00240141673353440, 
                       0.0240484422619185, 0.0407889947889709, 0.0329460704164021, 0.0131079241989148, 
                       0.00802463252306346, -0.00368723456237212, -0.0231357301749606, 
                       -0.026610152022363, -0.0159094648851324, -0.0188213074033152, 
                       -0.00970335155584632, -0.021235064079561, -0.0230462009765107, 
                       0.0133111894488440, 0.0148031665581324, -0.00603229151436025, 
                       3.46233919938086e-05, -0.0108024827876326, -0.000932596310983757, 
                       -0.0169929384781775, -0.00992393699591924, 0.0132764364040259, 
                       0.00205627881925053, 0.0258112340481444, 0.0281776480654252, 
                       0.0160149523840089, 0.00523373705029034, -0.0169637836588420, 
                       0.000623112029078655, -0.0141868755826321, -0.00321171121044831, 
                       0.0743950495497493, 0.0542289015672049, 0.0391968714568537, -0.0135541920624836, 
                       0.0153542242613047, -0.0094194197841322, -0.00547094071553647, 
                       -0.00239913220048973, 0.019879169496232, 0.0166006937251467, 
                       -0.00884289714417078, -0.0149238663033642, -0.0423907724837532, 
                       -0.0381010011953485, -0.00911201302570102, -0.0169016146575456, 
                       0.0105083286881965, 0.058273008019799, 0.0434196036832031, -0.0215951469910742, 
                       -0.00230126625611074, 0.00547170676692825, -0.00880871807226957, 
                       -0.0464745776581838, -0.0325544204388365, 0.0131149403889219, 
                       0.0033492489048061, 0.00754822245041555, -0.0097024181870838, 
                       -0.00965382308429241, -0.0100095947543497, -0.0183942240932586, 
                       0.0122552165022652, 0.0297247068675381, 0.0218407242637185, -0.0104993749706284, 
                       -0.0292879522068121, -0.0187554079279893, -0.0255819646811051, 
                       -0.0394480009717024, -0.0293045410816044, 0.0102755363948094, 
                       0.0230202055916129), timeindex = as.integer(c(13, 14, 22, 24, 
                                                                     25, 27, 30, 41, 89, 95, 4, 20, 24, 47, 49, 52, 61, 62, 69, 91, 
                                                                     3, 14, 32, 34, 45, 53, 54, 71, 83, 87, 3, 8, 18, 28, 41, 46, 
                                                                     50, 77, 92, 95, 12, 16, 20, 29, 39, 49, 57, 63, 75, 88, 9, 23, 
                                                                     26, 29, 34, 37, 39, 40, 45, 94, 26, 36, 39, 45, 50, 53, 62, 77, 
                                                                     91, 98, 11, 30, 31, 36, 46, 49, 57, 58, 90, 94, 2, 15, 30, 51, 
                                                                     64, 71, 77, 79, 87, 96, 30, 37, 38, 40, 57, 70, 76, 77, 89, 92, 
                                                                     3, 6, 20, 23, 41, 45, 53, 60, 89, 97, 11, 17, 21, 25, 30, 41, 
                                                                     63, 68, 84, 86, 4, 17, 22, 35, 52, 60, 61, 74, 88, 99, 11, 37, 
                                                                     43, 47, 56, 66, 72, 85, 95, 98, 14, 47, 50, 56, 59, 69, 79, 83, 
                                                                     87, 100, 10, 15, 25, 38, 48, 58, 61, 62, 79, 87, 10, 22, 37, 
                                                                     72, 79, 83, 84, 90, 96, 100, 28, 34, 38, 44, 47, 52, 54, 65, 
                                                                     73, 96, 13, 19, 26, 46, 48, 54, 55, 75, 85, 91, 1, 13, 16, 29, 
                                                                     37, 54, 58, 66, 72, 99, 12, 25, 34, 37, 54, 62, 70, 79, 84, 93, 
                                                                     34, 55, 59, 61, 73, 76, 77, 79, 86, 87, 1, 9, 12, 20, 25, 46, 
                                                                     56, 82, 86, 97, 6, 12, 43, 57, 60, 61, 69, 82, 84, 87, 1, 2, 
                                                                     6, 9, 16, 35, 49, 69, 85, 99, 18, 20, 21, 31, 34, 44, 46, 53, 
                                                                     70, 89, 13, 37, 40, 42, 44, 46, 55, 70, 72, 77, 5, 24, 31, 33, 
                                                                     42, 44, 52, 65, 66, 74, 13, 20, 21, 25, 36, 71, 78, 82, 89, 91, 
                                                                     3, 32, 35, 45, 66, 74, 83, 93, 94, 99, 5, 18, 23, 38, 50, 51, 
                                                                     72, 75, 91, 92, 5, 11, 14, 15, 19, 43, 49, 83, 87, 91, 19, 25, 
                                                                     30, 33, 34, 61, 81, 82, 87, 94, 1, 8, 16, 21, 25, 29, 69, 80, 
                                                                     97, 98, 10, 38, 58, 59, 65, 72, 73, 88, 90, 98, 14, 16, 26, 30, 
                                                                     31, 34, 59, 63, 74, 91, 9, 18, 37, 43, 65, 85, 91, 92, 96, 97, 
                                                                     4, 17, 18, 26, 32, 51, 57, 64, 81, 90, 14, 25, 30, 38, 51, 52, 
                                                                     54, 73, 82, 97, 6, 22, 34, 37, 40, 42, 55, 57, 69, 76, 16, 26, 
                                                                     35, 48, 55, 65, 67, 78, 81, 88, 7, 9, 12, 26, 36, 47, 48, 56, 
                                                                     84, 93, 16, 27, 28, 35, 40, 43, 55, 69, 83, 92, 20, 25, 27, 29, 
                                                                     38, 55, 72, 75, 82, 86, 11, 25, 49, 62, 67, 70, 79, 87, 91, 94, 
                                                                     4, 7, 10, 13, 16, 22, 29, 36, 58, 67, 1, 9, 28, 47, 52, 57, 72, 
                                                                     80, 90, 93, 3, 18, 19, 43, 60, 64, 68, 75, 81, 82, 8, 10, 11, 
                                                                     30, 32, 40, 59, 87, 99, 100, 1, 8, 18, 45, 68, 71, 77, 82, 88, 
                                                                     91, 3, 8, 10, 21, 29, 47, 61, 76, 83, 85, 18, 37, 44, 50, 52, 
                                                                     56, 84, 85, 88, 94, 3, 13, 14, 18, 23, 25, 58, 69, 87, 94, 1, 
                                                                     9, 12, 19, 64, 79, 81, 83, 85, 89, 13, 16, 20, 39, 46, 61, 64, 
                                                                     93, 94, 97, 17, 24, 26, 39, 61, 63, 68, 77, 80, 84, 17, 25, 46, 
                                                                     49, 56, 59, 60, 70, 75, 96, 22, 35, 55, 71, 77, 82, 85, 86, 93, 
                                                                     98, 21, 34, 36, 41, 47, 48, 54, 79, 80, 92, 7, 12, 14, 28, 36, 
                                                                     64, 74, 82, 88, 95, 15, 17, 32, 41, 61, 62, 68, 72, 80, 97, 10, 
                                                                     16, 27, 31, 38, 46, 57, 69, 72, 81, 19, 28, 34, 35, 47, 66, 67, 
                                                                     78, 88, 99, 20, 28, 47, 49, 50, 71, 88, 90, 93, 97, 18, 35, 53, 
                                                                     61, 65, 66, 72, 92, 93, 95, 5, 9, 15, 26, 40, 42, 62, 63, 75, 
                                                                     90, 3, 10, 12, 13, 32, 55, 65, 89, 92, 94, 4, 6, 16, 23, 45, 
                                                                     52, 78, 83, 96, 99, 20, 25, 37, 43, 46, 49, 54, 69, 70, 71, 3, 
                                                                     11, 19, 28, 34, 55, 61, 70, 74, 76, 3, 33, 36, 54, 55, 56, 60, 
                                                                     72, 88, 95, 1, 24, 27, 39, 49, 66, 72, 83, 95, 100, 19, 22, 30, 
                                                                     46, 51, 57, 60, 68, 72, 97, 5, 22, 38, 56, 60, 69, 79, 81, 91, 
                                                                     94, 4, 14, 50, 52, 72, 74, 92, 98, 99, 100, 26, 36, 37, 41, 46, 
                                                                     55, 60, 73, 74, 99, 4, 14, 19, 34, 43, 49, 69, 88, 93, 96, 11, 
                                                                     16, 35, 41, 53, 58, 62, 67, 69, 97, 9, 10, 19, 31, 41, 47, 50, 
                                                                     71, 89, 95, 1, 8, 14, 19, 22, 35, 43, 55, 77, 83, 7, 38, 45, 
                                                                     51, 66, 78, 79, 88, 93, 100, 1, 5, 64, 73, 74, 76, 86, 91, 92, 
                                                                     97, 22, 55, 62, 64, 66, 71, 75, 91, 92, 94, 8, 11, 25, 42, 58, 
                                                                     67, 72, 82, 84, 98, 3, 30, 33, 55, 67, 73, 76, 87, 92, 99, 7, 
                                                                     25, 29, 38, 45, 47, 66, 67, 81, 99, 2, 11, 18, 22, 25, 40, 51, 
                                                                     62, 90, 97, 23, 29, 40, 42, 48, 50, 60, 63, 89, 95, 1, 35, 37, 
                                                                     43, 53, 64, 69, 75, 77, 80, 3, 21, 30, 39, 41, 61, 65, 70, 81, 
                                                                     85, 16, 19, 23, 36, 47, 51, 57, 60, 73, 80, 4, 5, 8, 31, 66, 
                                                                     67, 71, 82, 89, 91, 9, 15, 23, 30, 38, 62, 65, 74, 84, 86, 1, 
                                                                     3, 27, 43, 45, 48, 53, 80, 85, 89, 9, 22, 25, 35, 36, 51, 73, 
                                                                     80, 85, 92, 26, 32, 39, 40, 48, 50, 56, 76, 78, 89, 2, 37, 43, 
                                                                     47, 65, 68, 71, 84, 92, 98, 4, 6, 7, 17, 22, 24, 54, 65, 86, 
                                                                     95, 6, 8, 45, 56, 58, 61, 62, 65, 88, 92, 4, 5, 33, 35, 39, 56, 
                                                                     57, 59, 83, 91)), curve = as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 
                                                                                                            1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
                                                                                                            3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
                                                                                                            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 
                                                                                                            8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 
                                                                                                            10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 
                                                                                                            11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 
                                                                                                            13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 
                                                                                                            14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 
                                                                                                            16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 
                                                                                                            18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 
                                                                                                            19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 
                                                                                                            21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 
                                                                                                            22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 
                                                                                                            24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 26, 
                                                                                                            26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 
                                                                                                            27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 
                                                                                                            29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 
                                                                                                            30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 
                                                                                                            32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 
                                                                                                            34, 34, 34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 
                                                                                                            35, 35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 37, 37, 37, 
                                                                                                            37, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 
                                                                                                            38, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, 40, 
                                                                                                            40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 42, 
                                                                                                            42, 42, 42, 42, 42, 42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 43, 
                                                                                                            43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, 45, 45, 
                                                                                                            45, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 46, 46, 46, 
                                                                                                            46, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 
                                                                                                            48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 
                                                                                                            50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 
                                                                                                            51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 
                                                                                                            53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 
                                                                                                            54, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 
                                                                                                            56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 
                                                                                                            58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, 59, 
                                                                                                            59, 59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 61, 61, 61, 
                                                                                                            61, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 62, 62, 62, 
                                                                                                            62, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 
                                                                                                            64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 66, 
                                                                                                            66, 66, 66, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 67, 67, 67, 
                                                                                                            67, 67, 67, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 69, 69, 69, 
                                                                                                            69, 69, 69, 69, 69, 69, 69, 70, 70, 70, 70, 70, 70, 70, 70, 70, 
                                                                                                            70, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 
                                                                                                            72, 72, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 74, 
                                                                                                            74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 
                                                                                                            75, 75, 75, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 77, 77, 77, 
                                                                                                            77, 77, 77, 77, 77, 77, 77, 78, 78, 78, 78, 78, 78, 78, 78, 78, 
                                                                                                            78, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 80, 80, 80, 80, 80, 
                                                                                                            80, 80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 82, 
                                                                                                            82, 82, 82, 82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, 83, 83, 
                                                                                                            83, 83, 83, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 85, 85, 85, 
                                                                                                            85, 85, 85, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 86, 86, 86, 
                                                                                                            86, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 88, 88, 88, 88, 88, 
                                                                                                            88, 88, 88, 88, 88, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 90, 
                                                                                                            90, 90, 90, 90, 90, 90, 90, 90, 90, 91, 91, 91, 91, 91, 91, 91, 
                                                                                                            91, 91, 91, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 93, 93, 93, 
                                                                                                            93, 93, 93, 93, 93, 93, 93, 94, 94, 94, 94, 94, 94, 94, 94, 94, 
                                                                                                            94, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 96, 96, 96, 96, 96, 
                                                                                                            96, 96, 96, 96, 96, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 98, 
                                                                                                            98, 98, 98, 98, 98, 98, 98, 98, 98, 99, 99, 99, 99, 99, 99, 99, 
                                                                                                            99, 99, 99, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)), 
                 gamma = structure(c(0.0218373389491933, -0.025668465747454, 
                                     0.152268529880326, 0.187235018524717, -0.0480609283955151, 
                                     0.0242005525149358, -0.0098284323142881, 0.0118000207663913, 
                                     0.188880028892866, -0.181724145628562, -0.0100510904485168, 
                                     -0.120663546655079, 0.0339437384177028, -0.227862892809654, 
                                     -0.192995563117666, 0.0673229970611408, 0.082071363284103, 
                                     0.154005999744787, 0.0340223295152608, -0.193589008262679, 
                                     0.00955236823720577, -0.105836745801777, -0.0209545973094202, 
                                     -0.0147644060467011, -0.0269140537149927, 0.0345489501454048, 
                                     -0.106995182431797, -0.0159237872625518, -0.0538136006461376, 
                                     0.151156389038849, -0.0327863736020349, 0.0857521480375477, 
                                     0.0374108624055816, 0.104655416579146, 0.0341585731802401, 
                                     -0.0392444060974585, -0.0924709654492949, 0.0423430600503454, 
                                     0.0984431552034132, -0.0779709150028069, 0.0116958827081058, 
                                     -0.0363224890368247, 0.138902448173654, -0.116299980311064, 
                                     -0.113223919169822, -0.0134021466896058, -0.0289130732462312, 
                                     -0.107550227573546, -0.0557010382352323, 0.171682938794451, 
                                     0.137662918080833, -0.0592549118253688, 0.0296671991011925, 
                                     0.176720901560753, -0.158052814284795, -0.0668991090835221, 
                                     0.106623006652369, -0.0655260101286987, 0.0286832620854542, 
                                     0.0899201572537778, 0.0382576881275488, 0.145210801741294, 
                                     0.0382345268858015, -0.101947736380216, 0.0524369369233153, 
                                     0.134413906866559, 0.0809765262960312, -0.156963199531822, 
                                     0.0262120195436506, -0.102117807868055, -0.0245758744288464, 
                                     0.12828492415034, -0.155620487019393, 0.244472487842511, 
                                     -0.166300269544407, -0.033777048463377, -0.058190350419559, 
                                     -0.248935512151962, 0.0388323155443128, 0.178888673154711, 
                                     0.0552506676688374, 0.103687113642831, 0.18128677799886, 
                                     0.083807877061334, 0.084925444255897, 0.0201332339563385, 
                                     -0.0244972782201977, 0.046902774077539, -0.086239588288116, 
                                     -0.0307706054317708, -0.145998852153208, 0.0930593846231661, 
                                     -0.120509839922316, 0.099953858833266, -0.0138947408743580, 
                                     -0.139471535134111, 0.0438922422408422, -0.0465987555897190, 
                                     -0.0155683948871406, -0.0228778196528604, -0.116657506602496, 
                                     -0.0398104500630008, 0.0115250330044457, -0.0305050798520915, 
                                     -0.0517944150268954, -0.0933915379252483, 0.0100309591127210, 
                                     -0.0131088194558165, 0.209322638526009, 0.149023619478461, 
                                     -0.094313956413595, 0.0494150360894934, -0.0396795209914649, 
                                     -0.117520547743341, 0.0336369924788886, -0.118606901183644, 
                                     -0.00561883804576216, 0.0238620529678781, 0.190263117118338, 
                                     0.0830704218916457, 0.133854394556927, -0.0194887780495371, 
                                     -0.00596702254784681, -0.0599653136393565, -0.128475565744945, 
                                     0.0799942986276874, -0.139310099532442, -0.113948335764725, 
                                     0.107970150184172, -0.096042677754784, -0.0762333061519201, 
                                     -0.0255103499222815, -0.0257153058706938, -0.0181719337700498, 
                                     0.0827880154412993, 0.0498791657500504, 0.0588073445171991, 
                                     0.132764017500600, 0.118398819992264, 0.0404664768068194, 
                                     -0.0726118934408867, 0.0759125892423778, -0.0399327067471508, 
                                     0.0154749052853122, -0.111851848059968, -0.111904802430703, 
                                     0.255230960694842, -0.0667125682241554, 0.149210234838749, 
                                     -0.0117126939302971, -0.128668811952505, -0.0275692259888875, 
                                     0.0552205725851597, -0.0159897863251010, -0.0525756613019199, 
                                     0.167268885695565, -0.0267832940036607, -0.00575313053954152, 
                                     -0.0193561893321833, 0.100195969024376, -0.0305997708194838, 
                                     -0.166981883341595, -0.0768301594131943, -0.05292541249248, 
                                     0.183701470458235, 0.117497227241963, -0.230216013226373, 
                                     0.091410832719235, -0.0748602415664677, 0.0465271259365904, 
                                     0.0694942834551174, -0.00628053077223831, 0.0874787353216655, 
                                     -0.0783758816912216, 0.114868858538408, 0.0647004550102879, 
                                     -0.063164583889129, -0.142956970126486, 0.0147040646470154, 
                                     0.0749171958239057, -0.233193871138992, -0.0130487973879127, 
                                     0.0397849830809272, -0.0109886931543078, 0.0175762064293732, 
                                     0.0960585192393703, 0.0546305673030516, -0.00955120491839235, 
                                     -0.0742449621200783, -0.217246362544907, 0.100018142758423, 
                                     -0.0232726733006344, 0.050618246627282, 0.0706261139107936, 
                                     0.0160075701594259, 0.0655394439317799, -0.0423108169075952, 
                                     -0.0579260078032337, -0.0313373986607216, 0.0561093751746339, 
                                     -0.0348154301243812, 0.00735317740875592, -0.16302576046613, 
                                     0.127223494411719, 0.0485240685596917, 0.0189124423389823, 
                                     -0.00931231160864541, -0.0915487261608617, -0.00578595968015532, 
                                     -0.0129709488353758, 0.104879972314105, 0.0086998760372691, 
                                     -0.0018700931472022, -0.123656763086121, -0.0215657902430287, 
                                     -0.0913519331172647, 0.0671908959704219, -0.0807818967868952, 
                                     0.0893879662008441, 0.000132911325147180, -0.0353250630284835, 
                                     -0.0300722208473015, -0.0362274857536355, 0.0193103773729370, 
                                     0.0622950831107744, -0.0193393792210944, -0.199329198402364, 
                                     -0.220451357616527, 0.177012462219202, 0.0250653602664899, 
                                     0.0277847684970977, -0.0301823436352857, -0.0429003413625956, 
                                     0.0100077860656836, -0.0381802396008532, -0.051683980313319, 
                                     0.0589366678798811, -0.136770623131164, -0.103676700739268, 
                                     -0.0772538266988803, 0.0188712310096444, 0.0209773447853462, 
                                     -0.0921808495726365, 0.0500141437887573, 0.05742653019906, 
                                     -0.148976824117844, -0.0151215523236331, -0.133547530414114, 
                                     -0.0214386380784278, -0.0369805983676153, -0.107947656050772, 
                                     0.146498569103698, 0.0961290108889362, 0.0145627856580497, 
                                     -0.0355762817778976, 0.0194710481118874, 0.0284562306193320, 
                                     0.094788129789424, -0.0236468255136769, -0.101823020825684, 
                                     -0.123664963361206, -0.021746500358619, -0.0621537686247511, 
                                     -0.145405942635031, 0.175108205270362, 0.160178614239986, 
                                     -0.041955122783418, -0.0571793106321322, -0.0918908419170996, 
                                     0.160830904248433, 0.0448585806051643, 0.0554327373178744, 
                                     -0.0401910634500659, 0.124638775453580, 0.0499130189764488, 
                                     0.0367663723913789, -0.0407325545833566, 0.00301538123055146, 
                                     -0.168162492651550, -0.0258216461460019, -0.0267187111515229, 
                                     0.00305197983091863, 0.0382989638773183, -0.0947774947452556, 
                                     0.0073913614222953, -0.137514843824587, -0.120869659584808, 
                                     -0.0522281370657635, 0.0414665488789786, -0.00523074513881978, 
                                     0.0735443814636765, -0.124791219178356, -0.0331857188881184, 
                                     0.0198694528557031, 0.0894215146191858, -0.0331159259138130, 
                                     -0.106812436470815, 0.0390553938387962, 0.0473292758825736, 
                                     0.166633232760316, -0.0762088074362277, -0.105367729977921, 
                                     -0.00617982432969952, -0.0529126500732911, -0.0489988883423116, 
                                     0.0734843383724082, 0.0407871492037317, -0.0098972048568866, 
                                     0.0316407438772775, 0.00923659273081106, -0.0250022782229929, 
                                     -0.0381233564953059, -0.0385780309439687, 0.116504792511739, 
                                     0.0737078330074053, -0.0117875275133418, -0.0694566188755314, 
                                     -0.00476556980145994, -0.0152309863301293, -0.0697714754235009, 
                                     -0.0460040045918093, 0.0456511247505918, 0.0647689171803875, 
                                     -0.0677724518451511, 0.146603414060990, -0.0322282588565605, 
                                     -0.0464914061845817, 0.143324065284417, 0.184753734958297, 
                                     -0.0288429113937392, -0.0634640949018798, -0.124922819437506, 
                                     0.141718113166246, -0.0787322065702254, -0.0167761715304091, 
                                     -0.230681855046164, -0.0138995110452554, 0.0906710086296502, 
                                     0.0605966575463631, -0.02332156369036, -0.0655659327858163, 
                                     0.166679110428269, 0.0045518026253642, 0.101649429624414, 
                                     0.0815829810294948, -0.10986433038531, 0.0312099469112287, 
                                     -0.0142310000642699, -0.00108943198281863, 0.100535409990239, 
                                     -0.04320421364032, 0.00271878230475911, 0.0507635477470265, 
                                     0.0522174398479201, 0.061680362433097, 0.0237598247439731, 
                                     0.239284218538152, 0.0332985953029462, -0.0696358069403424, 
                                     0.126233782877627, 0.146108946014264, -0.000764539089977187, 
                                     0.160894267191654, 0.155871061344269, 0.094005849908491, 
                                     -0.119471432910055, -0.176754680020080, -0.000913010832518756, 
                                     -0.206292872377721, 0.106352765028278, 0.0541706652764951, 
                                     0.108791993988535, -0.0355580891791312, -0.125152318594200, 
                                     0.0457955070620411, 0.0242184914808950, 0.108655845304490, 
                                     0.02727499544564, 0.172771536029315, -0.00626876478398636, 
                                     -0.0275641289613861, 0.0443881733930832, 0.0150724644030881, 
                                     -0.0101433540198459, -0.0157113915444035, -0.027864509763848, 
                                     -0.0304870643770814, -0.00174096998847598, 0.0044963079779712, 
                                     0.00410428626658667, 0.107452432806718, -0.113589803168452, 
                                     -0.0569230540285154, 0.0706050480798885, 0.00851142792290795, 
                                     -0.208326826928613, -0.0428524700717186, -0.0707921563685498, 
                                     0.212752379468486, -0.0306959627355285, -0.0942519952133903, 
                                     0.0730399122215053, -0.00322471294338434, -0.185958938909656, 
                                     -0.0982009064136909, 0.0316388649696724, -0.163052502691546, 
                                     0.113291921056650, 0.0918849346836849, 0.0613051837061325, 
                                     -0.137061443132983, 0.0141485630653158, -0.0935305674813998, 
                                     0.122559425759775, -0.0670021881951157, 0.00453089266391953, 
                                     0.0271795790243230, -0.0621333910731312, 0.00969690108586031, 
                                     0.0261292145835096, 0.136283270554558, 0.00600629092434467, 
                                     0.0108473032920974, 0.0614581323811488, -0.0475667908873933, 
                                     -0.0791245624382605, 0.0593680447102009, -0.151957778299106, 
                                     0.0401568504347897, -0.107313425022772, -0.103675388855501, 
                                     -0.0162114017241312, 0.0398447476178016, -0.136127904345563, 
                                     -0.0957381537709076, 0.183411241165530, 0.0784678139325162, 
                                     -0.141050786404220, 0.0239168272418172, -0.0287891259137622, 
                                     0.33866024703514, -0.180069348444646, -0.0827104192076405, 
                                     -0.0436007205646515, -0.040309295243175, 0.00840647065064834, 
                                     -0.0138556178346088, 0.0433647145968717, 0.0848789083014942, 
                                     0.0450894967892586, -0.195847179989961, -0.160273174931598, 
                                     -0.0431471332922125, 0.00118480780702261, -0.146086242711495, 
                                     -0.107676406915319, -0.0975046904707011, 0.0518401555463257, 
                                     -0.0191498030329502, 0.0457192901684552, -0.124488734470833, 
                                     0.0176665150148298, 0.0856011364511653, -0.0457893871610213, 
                                     0.184828210886861, -0.000945352569386167, 0.138343039569128, 
                                     -0.327405971079352, 0.140663398407534, -0.156525488955882, 
                                     0.0907521902822644, 0.155318319084749, -0.138688488470347, 
                                     -0.0649352475394055, 0.0781251253002645, 0.144122066650274, 
                                     0.0950066624663367, 0.166119817443515, 0.0658015285359967, 
                                     -0.0464497746638864, 0.177066636081813, -0.00207001953330559, 
                                     -0.0417859099756309, -0.195474422300730, 0.0541967350634747, 
                                     0.0206984542671804, 0.0295919705456860, -0.219623979316872, 
                                     0.112484494182752, 0.136994678295922, 0.00365346490680398, 
                                     0.0676753844058685, 0.112247916179481, -0.0375975771518596, 
                                     -0.0608229730949876, 0.084665242991219, 0.0784919408005408, 
                                     0.124145866790008, 0.118931563147425, -0.0731325262332848
                 ), .Dim = as.integer(c(100, 5)))), .Names = c("x", "timeindex", 
                                                               "curve", "gamma"))
