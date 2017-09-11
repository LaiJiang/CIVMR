
#' @title linear algebra decompositions for CIV. (internal function.)
#' @description This function implements linear algebra steps to acquire necessary matrices for CIV construction.
#' @param G: original instruments with dimension nXp.
#' @param X: phenotype of interest. dimension nXk.
#' @param Z: possible pleiotropic phenotypes which have been measured. dimension nXr.
#' @keywords Cholesky decomposition
#' @return A list of matrices which will be called by solv_pcc() and pcc_IV().
#' @examples
#' data(simulation)
#' LA_decomposition(simulation$G,simulation$X,simulation$Z)
#' @export


LA_decomposition <- function(G,X,Z){

  temp_cor <- cor(G)
  diag(temp_cor) <- 0
  if(max(abs(temp_cor))>0.99)print("redundant SNPs in G")

  M <- X %*% solve(t(X)%*%X) %*% t(X)

  temp <- t(G)%*%G
  e <- eigen(temp)
  eigen_value <- e$values
  #to enforce nonsingular eigenvalue
  eigen_value[which(e$values < 0)] <- 0
  temp <- e$vectors

  if(FALSE){
    J = Diagonal( x= eigen_value )
    Tempinv = solve(temp)
    temp %*% J %*% Tempinv
  }


  GG_square <- temp %*% diag(sqrt(eigen_value)) %*% t(temp)
  #GG_square is the square root of matrix G'G
  inv_GG_square <- solve(GG_square)


  #T <- tcrossprod(G,inv_GG_square)
  T <- inv_GG_square
  trans_T <- t(T)
  #T is the transformation matrix Tc = Gu

  #the problem is translated into max c'Ac with c'c=1 and B'c = 0
  A <- inv_GG_square %*% t(G) %*% M %*% G %*% inv_GG_square
  B <- inv_GG_square %*% t(G) %*% Z

  list(A=A,B=B,M=M, GG_square=GG_square, inv_GG_square=inv_GG_square,T=T, trans_T = trans_T)
}





#' @title Find a unique solution of CIV (internal use).
#' @description This function find a unique solutin to the constrained instrument problem given matrix A and B.
#' @param A: matrix given by LA_decomposition().
#' @param B: matrix given by LA_decomposition().
#' @keywords Cholesky decomposition
#' @return c: solution to the constrained maximization problem.
#' @return max_value: the maximized correlation value.
#' @examples
#' data(simulation)
#' #CIV linear algebra decomposition components
#' civ.deco <- LA_decomposition(simulation$G,simulation$X,simulation$Z)
#' #solve the CIV solution \eqn{u}
#' civ.c <- solve_pcc(civ.deco$A, civ.deco$B)
#' #plot the weight c
#' plot(civ.c$c)
#' @export
solve_pcc <- function(A,B){

  p <- nrow(B)
  k <- ncol(B)

  qrd <- qr(B)
  Q <- qr.Q(qrd,complete= TRUE)
  R <- qr.R(qrd,complete= TRUE)

  QTAQ <- crossprod(Q,crossprod(A, Q))

  A22 <- QTAQ[(k+1):(p),(k+1):(p)]

  eg <- eigen(A22)

  d <- eg$vectors[,1]


  c <- Q %*% c(rep(0,k),d)

  #the maximized value of correlation
  CTAC <- crossprod(c,crossprod(A,c))

  list(c=c,max_value=CTAC)
}



#' @title multiple orthogonal CIV solutions. (internal function)
#' @description This function find multiple CIV solutions that are orthogonal to each other. Only the first one achive the global maximum correlation.
#' @param A: matrix given by LA_decomposition().
#' @param B: matrix given by LA_decomposition().
#' @param G: original instruments.
#' @return u_max: the solution of u that would maximize the constrained correlation problem.
#' @examples
#' data(simulation)
#' #CIV linear algebra decomposition components
#' civ.deco <- LA_decomposition(simulation$G,simulation$X,simulation$Z)
#' #solve the CIV solution for c
#' civ.mult <- pcc_IV(civ.deco$A, civ.deco$B, simulation$G, civ.deco$inv_GG_square)
#' @export
pcc_IV <- function(A,B,G, inv_GG_square,no_IV=ncol(G)-ncol(B)){

  current_B <- B
  iv_list <- NULL
  cor_list <- NULL

  p <- ncol(G)
  k <- ncol(B)
  for(i in 1:(no_IV)){

    pcc_result <- solve_pcc(A, current_B)
    current_c <- pcc_result$c
    current_max_value <- pcc_result$max_value

    iv_list <- cbind(iv_list,current_c )
    cor_list <- c(cor_list,current_max_value)

    current_B <- cbind(current_B,current_c)




  }
  #Output: the IVs Gu,  correlation list of all IVs
  #all IVs u_all, the strongest IV u_max
  u_all <-  (inv_GG_square %*%  iv_list  )
  list( Gu = G %*% u_all ,
        cor_list = cor_list, u_all =  u_all,
        u_max = u_all[,1] )

}


#' @title Find a unique solution of CIV.
#' @description This function find the unique CIV solution.
#' @param MR.data: A data.frame() object containg G,X,Z,Y.
#' @keywords LA_decomposition, solve_pcc.
#' @return c: solution vector to the constrained maximization problem.
#' @return max_value: the maximized correlation value.
#' @return CIV: the new instrumentable variable for MR analysis: constrained instrumental variable.

#' @examples
#' data(simulation)
#' civ.fit <- CIV(simulation)
#' #plot the weight c
#' plot(civ.fit$c)
#' @export
CIV <- function(MR.data){

  civ.deco <- LA_decomposition(MR.data$G,MR.data$X,MR.data$Z)

  solution <- solve_pcc(civ.deco$A, civ.deco$B)

  c <- solution$c

  new.G <- (MR.data$G) %*% c

  list(c=c,max_value=solution$max_value, CIV=new.G)
}



#' @title cross-validated CIV.
#' @description This function produce a Constrained Instrumental Variable with cross-validation.
#'              Specifically, for a predefined fold the CIV is calculated using
#'               all samples except this fold, then the CIV solution is applied to the samples  in this fold
#'               to obtain corresponding CIV. In this way the correlation between samples
#'               are expected to be reduced.
#' @param MR.data: a data frame containing G,X,Z,Y.
#' @param n_folds: number of folds for cross-validation.
#' @return weights: A matrix with dimension \eqn{n_folds * p}. Each row is a CIV solution \eqn{c} from a specific fold.
#' @return civ.IV: cross-validated CIV instrument \eqn{G^{*}=Gc}.
#' @return beta_est: causal effect estimation of X on Y using CIV instrument civ.IV
#' @examples
#' data(simulation)
#' cv.civ <- cv_CIV(simulation)
#' #strong correlation between CIV solutions from different folds.
#' cor(t(cv.civ$weights))
#' @export
cv_CIV <- function(MR.data, n_folds = 10 ){

  G <- MR.data$G
  X <- MR.data$X
  Z <- MR.data$Z
  Y <- MR.data$Y

  if( length(ncol(X))==0 ) X <- matrix(X,ncol=1)

  n <- length(Y)
  p <- ncol(G)

  #allele_score is used to save the calculated allele scores
  civ_score <- X
  #weights is a matrix where each row contains weight from each folds
  weights <- NULL
  flds <- createFolds(c(1:n), k = n_folds, list = TRUE, returnTrain = FALSE)


  for(cv_folds in 1:n_folds){

    apply_index <- flds[[cv_folds]]
    train_index <- c(1:n)[-apply_index]

    if(ncol(X)==1){
      X_train <- as.matrix(X[train_index,1],ncol=1)
      X_apply <- as.matrix(X[apply_index,1],ncol=1)
    }
    else{
      X_train <- X[train_index,]
      X_apply <- X[apply_index,]
    }

    #the split of training and test data
    G_train <- G[train_index,]
    G_apply  <- G[apply_index,]

    Z_train <- Z[train_index,]
    Z_apply <- Z[apply_index,]

    #now run the CIV training and find CIV weight
    deco <- LA_decomposition(G_train, X_train, Z_train)

    pcc_iv <- pcc_IV( deco$A, deco$B, G_train,
                      deco$inv_GG_square,no_IV = 1)

    #CIV training weight
    civ.u <- pcc_iv$u_max



    if(ncol(X)==1){
      if(ncol(G)>1){
        civ_score_apply <- G_apply%*%civ.u}
      else{
        civ_score_apply <- G_apply%*%civ.u
      }
      civ_score[apply_index] <- civ_score_apply
    }




    weights <- rbind(weights,civ.u)
  }


  #now we use the cross-validated CIV score to do causal effect estimation
  MR.data <- MR.data
  MR.data$G <- civ_score

  est_TSLS <- TSLS_IV(MR.data)
  beta_est <- est_TSLS$coef[-1]


  list(beta_est = beta_est, weights = weights, civ.IV = civ_score)
}


#' @title bootstrapped CIV (recommended).
#' @description This function generate a bootstrapped CIV w/wo correction.
#'              Specifically, for a bootstrap sample we can generate civ solution u. The boostrap corrected
#'              solution u is obtained as the global solution u - ( bootrapped average u - global u).
#' @param MR.data: a data frame containing G,X,Z,Y.
#' @param n_boots: number of bootstrap samples.
#' @return boots.u: bootstrapped CIV solution u (without correction).
#' @return boots.cor.u: bootstrap corrected solution of u. (suggested)
#' @examples
#' data(simulation)
#' boot.civ <- boot_CIV(simulation)
#' #plot the bootstrap corrected solution u.
#' plot(boot.civ$boots.cor.u)
#' @export


boot_CIV <- function(MR.data, n_boots = 10 ){

  X <- MR.data$X
  G <- MR.data$G
  Z <- MR.data$Z
  Y <- MR.data$Y

  if( length(ncol(X))==0 ) X <- matrix(X,ncol=1)

  n <- length(Y)
  p <- ncol(G)


  civ.boot.mat <- NULL

  for(i_boot in 1:n_boots){

    obs_id <- sample(1:n, size=n, replace=TRUE)

    boot_X <- X[obs_id,]
    boot_Z <- Z[obs_id,]
    boot_G <- G[obs_id,]
    boot_Y <- Y[obs_id]

    #to avoid singularity
    cor_g <- cor(boot_G)
    #cor_z <- cor(boot_Z)
    diag(cor_g) <- 0

    #we keep re-sampling if the previous sampled Gs are highly correlated
    while( (is.na(max(cor_g))==1)||(max(abs(cor_g))>0.999) ){

      obs_id <- sample(1:n, size=n, replace=TRUE)

      boot_X <- X[obs_id,]
      boot_Z <- Z[obs_id,]
      boot_G <- G[obs_id,]
      boot_Y <- Y[obs_id]

      cor_g <- cor(boot_G)
      #cor_z <- cor(boot_Z)
      diag(cor_g) <- 0
    }

    boot_MR.data <- data.frame(Y = boot_Y,
                               n=n, p=p)
    boot_MR.data$G <- boot_G
    boot_MR.data$Z <- boot_Z
    boot_MR.data$X <- boot_X


    #
    deco <- LA_decomposition(boot_MR.data$G,boot_MR.data$X,boot_MR.data$Z)
    pcc_iv <- pcc_IV(deco$A,deco$B,boot_MR.data$G,
                     deco$inv_GG_square,no_IV = 1)
    #weight
    civ.u <- pcc_iv$u_max

    if(is.complex(civ.u)==0){civ.boot.mat <- cbind(civ.boot.mat,civ.u)}

  }


  deco <- LA_decomposition(MR.data$G,MR.data$X,MR.data$Z)
  pcc_iv <- pcc_IV(deco$A,deco$B,MR.data$G,
                   deco$inv_GG_square,no_IV = 1)
  #weight
  civ.global.u <- pcc_iv$u_max

  boots.u <- apply(civ.boot.mat,1,mean)

  boots.cor.u <- 2*civ.global.u - boots.u

  if(is.complex(boots.cor.u)==0){boots.cor.u <- boots.u}

  list(boots.u = boots.u, boots.cor.u = boots.cor.u)
}



#' @title SNP pre-processing.
#' @description This function remove highly correlated SNPs. It also
#' calculate MAF for each snp and delete rare snps with low MAF (e.g. 0.01).

#' @param snp_matrix: SNP matrix with dimension n * p.
#' @param crit_high_cor: criteria to choose highly correlated SNPs.
#' @param maf_crit: criteria to choose rare SNPs.
#' @return sel_snp: the new dosage matrix of selected SNPs.
#' @return id_snp: the ids (columns) of selected SNPs in the original SNP matrix.
#' @examples
#' data(simulation)
#' snp.rdc <- SNP_reduction(simulation$G)
#' @export


SNP_reduction <- function(snp_matrix,crit_high_cor=0.8,maf_crit =0.01){

  #remove high correlated SNPs

  #delete all redundant SNPs and repeating SNPs.
  sel_snp <- snp_matrix[,1]
  id_snp <- 1
  for(j in 2:ncol(snp_matrix)){
    if( (max(abs(cor(sel_snp,snp_matrix[,j]) )) <crit_high_cor) && (length(unique(snp_matrix[,j]))>1)  )
    {sel_snp <- cbind(sel_snp, snp_matrix[,j])
    id_snp <- c(id_snp,j)
    }
  }
  ########################
  if( length(ncol(sel_snp))==0){sel_snp <- matrix(sel_snp,ncol=1)}
  ########################
  #calculate MAF for each snp and delete rare snps with 0.01 MAF
  N <- nrow(sel_snp)
  mat_freq <- NULL
  for(j in 1:ncol(sel_snp)){
    temp <-  c( length(which(sel_snp[,j]==0))/N, length(which(sel_snp[,j]==1))/N,
                length(which(sel_snp[,j]==2))/N )
    mat_freq <- rbind(mat_freq,temp)
  }
  MAF <- apply(mat_freq[,c(2,3)],1,min)
  non.rare.snp <- which(MAF>maf_crit)
  sel_snp <- sel_snp[,non.rare.snp]
  id_snp <- id_snp[non.rare.snp]

  list(sel_snp =sel_snp, id_snp=id_snp)

}



#' @title Instrumental variable reduction.
#' @description This function remove highly correlated IVs. An upgraded function SNP_reduction() is suggested.
#' @param snp_matrix: IV matrix with dimension n * p.
#' @param crit_high_cor: criteria to choose highly correlated SNPs. default is 0.8 correlation.
#' @return sel_snp: the selected IVs.
#' @return id_snp: the ids (columns) of selected IVs in the original IV matrix.
#' @examples
#' data(simulation)
#' snp.rdc <- IV_reduction(simulation$G)
#' @export

IV_reduction <- function(snp_matrix,crit_high_cor=0.8){


  #delete all redundant SNPs and repeating SNPs.
  sel_snp <- snp_matrix[,1]
  id_snp <- 1
  for(j in 2:ncol(snp_matrix)){
    if( (max(abs(cor(sel_snp,snp_matrix[,j]) )) <crit_high_cor) && (length(unique(snp_matrix[,j]))>1)  )
    {sel_snp <- cbind(sel_snp, snp_matrix[,j])
    id_snp <- c(id_snp,j)
    }
  }
  list(sel_snp =sel_snp, id_snp=id_snp)

}

#' @title CIV_smooth solution given \eqn{\lambda}. (Internal function)
#' @description This function finds a CIV_smooth solution of u given a value of \eqn{\lambda}. This function is mostly for internal use. smooth_CIV() is suggested for users to obtain optimal solutions of CIV_smooth.
#' @param initial: the initial point of u for updating. The CIV solution will be used as the initial point if no choice is made.
#' @param G: SNP matrix with dimension \eqn{n \times p}.
#' @param X: phenotype of interest.
#' @param Z: pleiotropic phenotype Z.
#' @param GTG: \eqn{G`G}
#' @param GTMG: \eqn{G`X(X`X)^{-1}X`G}.
#' @param ZTG: \eqn{Z`G}
#' @param GTZ: \eqn{G`Z}
#' @param ZTG_ginv: general inverse of \eqn{Z`G} (ginv(\eqn{Z`G})).
#' @param null_space: null space of matrices G`Z (null(G`Z)).
#' @param lambda: a given value (must be specified) for regularization parameter \eqn{\lambda}.
#' @param accuracy_par: the accuracy threshold parameter to determine if the algorithm converged to a local maximum. Default is 1e-10.
#' @param last_conv_iters: the maximum iterations to run. Default is 2000.
#' @param ......: default values for other tuning parameters.
#' @return mat_u: the trace of all updated iterations of u.
#' @return opt_solution: the final solution of u.
#' @return value_list: the iteration values of target function (penalized correlation).
#' @return unstrained_val_list: the iteration values of correlation between X and Gu.
#' @return dev_list: the iteration values of deviance between updated vector of u.
#' @return n_iters_stage: the number of iterations before finishing updating. If this value < last_conv_iters, then the algorithm stopped at a solution of u without using up its updating quota.
#' @return sigma_stage: the updating values of \eqn{\sigma} that are used in each iteration.
#' @return stepsize_list: the updating values of stepsize that are used in each iteration.
#' @examples
#' data(simulation)
#' G <- simulation$G
#' X <- simulation$X
#' Z <- simulation$Z
#' GTG <- crossprod(G,G)
#' M <- tcrossprod ( tcrossprod ( X , solve(crossprod(X,X) ) ), X )
#' GTMG <- crossprod(G, crossprod(M,G))
#' ZTG <- crossprod(Z,G)
#' GTZ <- crossprod(G,Z)
#' null_space <- Null( GTZ)
#' ZTG_ginv <- ginv(ZTG)
#' lambda <- 1
#' smooth.lambda1 <- smooth_L0_lambda(null_space = null_space, G = G, X = X, GTG = GTG, lambda = lambda,
#' GTMG = GTMG, ZTG = ZTG, GTZ = GTZ, ZTG_ginv = ZTG_ginv )
#' plot(smooth.lambda1$opt_solution)  #plot the final solution u
#' @export


smooth_L0_lambda <- function(initial = NULL, null_space, G,X,GTG,
                             lambda, sigma_min = 0.01, sigma_up =0.5,
                             stepsize=0.1, conv_iters = 5, stepsize_last = 0.0001,
                             last_conv_iters = 2000,
                             GTMG, ZTG, GTZ, ZTG_ginv, accuracy_par=1e-10 ){

  p <- ncol(G)
  n <- nrow(G)
  #initialization
  #the initial u set as an "average" vector in null space
  #there is no randomness in this algorithm so for different lambda
  #we always start from the same initial point
  if(length(initial)==0){
    #initial <- null_space %*% c( rnorm(ncol(null_space) ,mean=0,sd= 1) )

    #if we have multiple columns in X, we use sparse cca to find an initial guess
    if(length(ncol(X))>1){
      sparseCCA <- CCA(G,X)
      initial <- sparseCCA$u
    }

    else{
      if(p<n){
        deco <- LA_decomposition(G,X,Z)
        pcc_iv <- pcc_IV(deco$A,deco$B,G,
                         deco$inv_GG_square,no_IV = 1)
        initial  <- as.vector(pcc_iv$u_max)
      }

      else{
        #initial choice: LASSO or simple LR?
        #cvres<-cv.lars(G,X,K=10,type='lasso',max.steps=20)
        cvres <- cv.glmnet(G,X)
        #to avoid which.min(cvres$cv)==1 or 0.
        sAtBest<-cvres$lambda.1se
        res <- glmnet(G,X,family="gaussian",lambda=sAtBest)
        initial <- as.vector(res$beta)
      }
    }

  }
  #project back to constrained set
  initial <- initial - tcrossprod(ZTG_ginv, crossprod(initial, GTZ) )
  L2scale <- 1/ as.numeric(sqrt(crossprod(initial,crossprod(GTG,initial ) )) )
  initial <- initial*L2scale

  #decreasing order of sigma
  sigma <- 2*max(sum(abs(initial)))

  #print(sigma)

  u <- initial

  #initialize the stepsize
  stepsize_iters <- stepsize

  stepsize_list <- NULL
  target_value <-  lambda*p-lambda*sum(exp(-u^2/2/sigma^2)) -  crossprod(u,crossprod(GTMG, u))

  direction <- (lambda*u/(sigma^2))*exp(-u^2/2/(sigma^2) ) - 2*as.numeric(crossprod(u,GTMG))
  u1 <- u - stepsize_iters*direction
  temp_u <- u1/as.numeric(sqrt(crossprod(u1,crossprod(GTG,u1 ) )) )
  target_cor <- lambda*p-lambda*sum(exp(-temp_u^2/2/sigma^2)) -  crossprod(temp_u,crossprod(GTMG, temp_u))

  dev_iter <- NULL
  n_iters_stage <- NULL
  sigma_stage <- NULL
  u_mat <- NULL
  u_mat <- cbind(u_mat, u)

  while(sigma > sigma_min){

    iters <- 0
    temp_dev <- 10000

    while( (iters < conv_iters)&&(temp_dev>accuracy_par)  ){

      prev_u <- u
      #the gradient direction
      direction <- (lambda*u/(sigma^2))*exp(-u^2/2/(sigma^2) ) -
        2*as.numeric(crossprod(u,GTMG))

      u1 <- u - stepsize_iters*direction
      #back to Au=0 null space
      u <- u1 - tcrossprod(ZTG_ginv, crossprod(u1, GTZ) )

      #back to L2 constraint set
      L2scale <- 1/ as.numeric(sqrt(crossprod(u,crossprod(GTG,u ) )) )
      u <- u*L2scale



      #we record the projected u.
      u_mat <- cbind(u_mat, u)
      temp_dev <- sqrt(sum((abs( abs(u) - abs(prev_u) ))^2))/p

      dev_iter <- c( dev_iter, temp_dev )
      #check the target function value
      target_value <-  c( target_value, lambda*p-lambda*sum(exp(-u^2/2/sigma^2)) -  crossprod(u,crossprod(GTMG, u)) )

      temp_u <- u1/as.numeric(sqrt(crossprod(u1,crossprod(GTG,u1 ) )) )
      temp_cor <- lambda*p-lambda*sum(exp(-temp_u^2/2/sigma^2)) -  crossprod(temp_u,crossprod(GTMG, temp_u))
      target_cor <- c(target_cor, temp_cor  )

      #Bold driver strategy to get around the local minimum trap: we reduce the stepsize by 50% if the target value inreased
      #bigger than 10^(-10). We increase the stepsize by 5% if the target value (error) is decreased.
      temp_count <- length(target_value)

      if(target_value[temp_count] - target_value[temp_count-1] > accuracy_par ){stepsize_iters <- stepsize_iters *0.5 }

      if(target_value[temp_count] - target_value[temp_count-1] < 0 ){stepsize_iters <- stepsize_iters *1.05 }

      iters <- iters + 1

      stepsize_list <- c(stepsize_list,stepsize_iters)

    }
    n_iters_stage <- c(n_iters_stage, iters)
    sigma_stage <- c(sigma_stage, sigma)

    sigma <- sigma * sigma_up
  }

  #####################################################################
  #after converge find a better u solution
  sigma <- sigma/sigma_up
  stepsize_last <- stepsize*stepsize_last
  stepsize_iters <- stepsize_last

  temp_dev <- 100000
  iters <- 0

  #if it doesn't converge automatically we let it run for a minimum iterations: last_conv_iters
  while( (temp_dev >accuracy_par )&&( iters <  last_conv_iters )  ){

    prev_u <- u

    #the gradient direction
    direction <- (lambda*u/(sigma^2))*exp(-u^2/2/(sigma^2) ) -
      2*as.numeric(crossprod(u,GTMG))

    u1 <- u - stepsize_last*direction


    #back to Au=0 null space
    u <- u1 - tcrossprod(ZTG_ginv, crossprod(u1, GTZ) )

    #back to L2 constraint set
    L2scale <- 1/ as.numeric(sqrt(crossprod(u,crossprod(GTG,u ) )) )
    u <- u*L2scale

    #we replace temp_cor with the target value of u1
    temp_u <- u1/as.numeric(sqrt(crossprod(u1,crossprod(GTG,u1 ) )) )
    temp_cor <- lambda*p -lambda*sum(exp(-temp_u^2/2/sigma^2)) - crossprod(temp_u,crossprod(GTMG, temp_u))
    #temp_cor <- crossprod(u,crossprod(GTMG, u))

    #check the target function value
    target_value <-  c( target_value, lambda*p -lambda*sum(exp(-u^2/2/sigma^2)) - crossprod(u,crossprod(GTMG, u))  )

    u_mat <- cbind(u_mat, u)
    target_cor <- c(target_cor,  temp_cor )

    #temp_dev <-   sum( abs(prev_u) - abs(u)  )
    temp_dev <- sqrt(sum(( abs( abs(u) - abs(prev_u)))^2))/p

    dev_iter <- c( dev_iter, temp_dev )

    #the bolde driver strategy
    temp_count <- length(target_value)

    if(target_value[temp_count] - target_value[temp_count-1] > accuracy_par ){stepsize_iters <- stepsize_iters *0.5 }

    if(target_value[temp_count] - target_value[temp_count-1] < 0 ){stepsize_iters <- stepsize_iters *1.05 }

    iters <- iters + 1

    stepsize_list <- c(stepsize_list,stepsize_iters)


  }
  n_iters_stage <- c(n_iters_stage, iters )

  sigma_stage <- c(sigma_stage, sigma)


  #how to specify the optimal solution: two ways: the final converged solution. or the local minimum solution
  #in the last stage.

  #opt_solution <- u_mat[,ncol(u_mat)]

  #we access the target values of the u from the last stage
  last_stage_value <-   tail(target_value, n=  sum(  tail( n_iters_stage, 2  )  )   )
  opt_solution <- u_mat[,  (  length(target_value)-  sum(tail(n_iters_stage, 2))) + which.min(last_stage_value)]

  #output
  list(mat_u = u_mat, value_list = target_value, unstrained_val_list = target_cor,
       dev_list = dev_iter, n_iters_stage = n_iters_stage, sigma_stage = sigma_stage,
       stepsize_list = stepsize_list, opt_solution = opt_solution )

}

#' @title CIV_smooth solution with cross-validation.(recommended)
#' @description This function first find the optimal value of \eqn{\lambda} according
#' to projected prediction error with cross-validation. Then for a given \eqn{\lambda} value multiple intial
#' points are used to explore potentially multiple modes.
#' @param initial: the initial value for updating u.
#' @param G: SNP matrix with dimension \eqn{n \times p}.
#' @param X: phenotype of interest.
#' @param Z: pleiotropic phenotype Z.
#' @param Y: the disease outcome Y.
#' @param lambda_list: a list of values for regularization parameter lambda. A default list will be chosen if not provided.
#' @param k_folds: number of folds for cross-validation (to find optimum \eqn{\lambda}). default = 10.
#' @param n_IV: the number of initial points chosen to explore potential multiple modes. The converged
#' solutions will be screened to delete redundant solutions. So the final solutions will be less or equal to n_IV. default = 100.
#' @param sigma_min: the minimum value of \eqn{\sigma} (corresponding to the closeast approximation of \eqn{L_0} penalty). default = 0.01.
#' @param sigma_up: the moving down multiplier. \eqn{\sigma_{j+1} = sigma_up \times \sigma_{j}}. default = 0.5.
#' @param stepsize: the stepsize to move solution u. default = 0.1.
#' @param conv_iters: the maximum steps to allow updating when a converged solution is found. default =5.
#' @param stepsize_last: When a converged solution is found with stepsize, we update this solution with a smaller stepsize to achive a more precise
#' local maximum solution. default = 0.0001.
#' @param last_conv_iters: the maximum iterations to run in the stage of ``refining" optimum solution. default = 2000.
#' @param ......: default values for other tuning parameters.
#' @return opt_lambda: the chosen optimum value of \eqn{\lambda} corresponding to the minimum projected prediction error (see paper).
#' @return IV_mat: the final matrix of CIV instruments with respect to the opt_lambda. Each column is a new instrument.
#' @return u_mat: the final CIV solutions of u with respect to the opt_lambda. Each column is a converged solution.
#' @return G_pred_error_list: the projected prediction error according to the list values of \eqn{\lambda}.
#' @return Pred_error_list: the prediction error according to the list values of \eqn{\lambda}.
#' @examples
#' data(simulation)
#' G <- simulation$G
#' X <- simulation$X
#' Z <- simulation$Z
#' Y <- simulation$Y
#' smooth.opt <- smooth_CIV( G,X,Z,Y, k_folds = 10)
#' plot(smooth.opt$u_mat[,1])  #plot a solution u.
#' @export
smooth_CIV <- function(G,X,Z,Y, lambda_list = NULL, k_folds =10,
                          sigma_min = 0.01, sigma_up =0.5,
                          stepsize=0.1, conv_iters = 5,
                          stepsize_last = 0.0001,
                          last_conv_iters = 2000,
                          method_lambda ="er",
                          n_IV = 100){


  n <- nrow(G)
  p <- ncol(G)
  GTG <-  crossprod(G,G)
  M <- tcrossprod ( tcrossprod ( X , solve(crossprod(X,X) ) ), X )
  GTMG <- crossprod(G, crossprod(M,G))
  ZTG <- crossprod(Z,G)
  GTZ <- crossprod(G,Z)
  null_space <- Null( GTZ)
  ZTG_ginv <- ginv(ZTG)

  if(length(lambda_list)==0){lambda_list <- c(seq(from=0.01, to=0.09, by = 0.01),
                                              seq(from=0.1,to=1, by=0.1))}

  #for lambda_list we obtain the R2 list and MSE list
  Pred_error_list <- NULL
  G_pred_error_list <- NULL

  #############################
  #initialization
  #if we have multiple columns in X, we use sparse cca to find an initial guess

  #####if only 1 variable in X.
  #we use cross-validated LASSO to find initial estimation of u

  if(p<n){
    #pcc_iv <- pcc_IV(A,B,G, inv_GG_square)
    #initial  <- as.vector(pcc_iv$u_max)

    #we use the adjusted G
    initial <-  rep(1,p)
  }

  else{#may LASSO maybe LR?
    cvres <- cv.glmnet(G,X)
    #to avoid which.min(cvres$cv)==1 or 0.
    sAtBest<-cvres$lambda.1se
    res <- glmnet(G,X,family="gaussian",lambda=sAtBest)
    initial <- as.vector(res$beta)
  }



  #project back to constrained set
  initial <- initial - tcrossprod(ZTG_ginv, crossprod(initial, GTZ) )
  L2scale <- 1/ as.numeric(sqrt(crossprod(initial,crossprod(GTG,initial ) )) )
  initial <- initial*L2scale

  #preserve the global optimal initial value for the final smooth fitting
  global_initial <- initial

  pb <- txtProgressBar(min=0, max=length(lambda_list), style=3)
  finish <- 0

  flds <- createFolds(c(1:n), k = k_folds, list = TRUE, returnTrain = FALSE)



  for(lambda in lambda_list){

    #for each value of lambda we use crossed-validated dataset to
    #find optimal u with penalized correlation or final prediction error
    #or
    #or one standard error rule (see CrossValidationStandardDeviation.pdf)

    Pred_error <- 0
    G_pred_error <- 0

    for(cv_folds in 1:k_folds){

      test_index <- flds[[cv_folds]]
      #print(train_index)
      #print("\n")
      #train_index <- sample(c(1:n),size=round(n/k_folds),replace=FALSE)
      train_index <- c(1:n)[-test_index]


      if(length(ncol(X))==0){X_train <- as.matrix(X[train_index],ncol=1)
      X_test <- as.matrix(X[test_index],ncol=1)
      }
      else{
        if((ncol(X))==1){X_train <- as.matrix(X[train_index,],ncol=1)
        X_test <- X[test_index,1]
        }

        else{ X_train <- X[train_index,]
        X_test <- X[test_index,]
        }
      }

      Y_test  <- Y[test_index]
      G_train <- G[train_index,]
      G_test  <- G[test_index,]
      Z_train <- Z[train_index,]
      GTG_train <- crossprod(G_train,G_train)
      GTZ_train <- crossprod(G_train,Z_train)
      null_space_train <- Null( GTZ_train)
      M_train <- tcrossprod ( tcrossprod ( X_train , solve(crossprod(X_train,X_train) ) ), X_train )
      GTMG_train <- crossprod(G_train, crossprod(M_train,G_train))
      ZTG_train <- crossprod(Z_train,G_train)
      ZTG_ginv_train <- ginv(ZTG_train)

      #if p_train<n_train then CIV initial point (on global data) is used.
      #if p_train>n_train then in the inner-loop either CCA or CIV used on local data
      #(depending on the dimension of G_train)
      #if(ncol(G_train)>nrow(G_train)){initial <- NULL}

      smooth <- smooth_L0_lambda(initial = initial, null_space = null_space_train,
                                 G = G_train, X = X_train, GTG=GTG_train,
                                 lambda=lambda, sigma_min = sigma_min,
                                 sigma_up =sigma_up,
                                 stepsize=stepsize, conv_iters = conv_iters,
                                 stepsize_last = stepsize_last,
                                 last_conv_iters = last_conv_iters,
                                 GTMG = GTMG_train, ZTG=ZTG_train, GTZ = GTZ_train,
                                 ZTG_ginv=ZTG_ginv_train)

      #we keep the solution
      u <- smooth$opt_solution

      IV_test <- G_test%*%u

      #now we caluclated penalized correlation and prediction error
      MR.data <- data.frame(Y_test)
      MR.data$G <- IV_test
      MR.data$X <- X_test
      MR.data$Y <- Y_test
      cv_tsls_fit <- TSLS_IV(MR.data)
      Pred_error_fld <- Y_test - cbind(1,X_test) %*% (cv_tsls_fit$coef)

      PG <- tcrossprod(IV_test,IV_test)/sum((IV_test)^2)
      G_pred_error_fld <- PG %*% Pred_error_fld

      #now add the prediction error in this fold to the total error
      Pred_error <- Pred_error + sum(Pred_error_fld^2)
      G_pred_error <- G_pred_error + sum(G_pred_error_fld^2)

    }


    Pred_error_list <- c(Pred_error_list, Pred_error )
    G_pred_error_list <- c(G_pred_error_list, G_pred_error)

    finish <- finish +1
    setTxtProgressBar(pb, finish)
  }
  cat("\n")

  #optimal lambda with respect to penalized correlation or G_pred_error
  if(method_lambda=="er"){opt_lambda <- lambda_list[which.min(G_pred_error_list)]}

  ###############################
  #find the optimal IV using  the global optimal initial values. this steps is to ensure we have at least 1 solution of u.
  smooth <- smooth_L0_lambda(initial = global_initial, null_space = null_space,
                             G=G, X=X, GTG=GTG,
                             lambda=opt_lambda, sigma_min = sigma_min,
                             sigma_up =sigma_up,
                             stepsize=stepsize, conv_iters = conv_iters,
                             stepsize_last = stepsize_last,
                             last_conv_iters = last_conv_iters,
                             GTMG = GTMG, ZTG=ZTG, GTZ=GTZ, ZTG_ginv=ZTG_ginv)

  u <- smooth$opt_solution
  IV_mat <- G %*% u
  u_mat <- u


  pb <- txtProgressBar(min=0, max=(n_IV-1), style=3)
  finish <- 0
  #find multiple optimal IV using  different random initial values.
  for(temp in 1:(n_IV-1)){
    smooth <- smooth_L0_lambda(initial = rnorm(p), null_space = null_space,
                               G=G, X=X, GTG=GTG,
                               lambda=opt_lambda, sigma_min = sigma_min,
                               sigma_up =sigma_up,
                               stepsize=stepsize, conv_iters = conv_iters,
                               stepsize_last = stepsize_last,
                               last_conv_iters = last_conv_iters,
                               GTMG = GTMG, ZTG=ZTG, GTZ=GTZ, ZTG_ginv=ZTG_ginv)

    u <- smooth$opt_solution
    IV <- G %*% u
    IV_mat <- cbind(IV_mat,IV)
    u_mat <- cbind(u_mat,u)

    finish <- finish +1
    setTxtProgressBar(pb, finish)

  }
  cat("\n")

  ###############################################
  #some of these u are local minimum and myabe we can let it converge to
  #some points.
  #we run it again and see if there is any difference.
  #we first run it with fixed lambda.

  #
  TS_IV_mat <- NULL
  TS_u_mat  <- NULL

  pb <- txtProgressBar(min=0, max=length(lambda_list), style=3)
  finish <- 0

  for(temp in 1:(ncol(u_mat))){
    smooth <- smooth_L0_lambda(initial = u_mat[,temp], null_space = null_space,
                               G=G, X=X, GTG=GTG,
                               lambda=opt_lambda, sigma_min = sigma_min,
                               sigma_up =sigma_up,
                               stepsize=stepsize, conv_iters = conv_iters,
                               stepsize_last = stepsize_last,
                               last_conv_iters = last_conv_iters,
                               GTMG = GTMG, ZTG=ZTG, GTZ=GTZ, ZTG_ginv=ZTG_ginv)

    u <- smooth$opt_solution
    IV <- G %*% u
    TS_IV_mat <- cbind(TS_IV_mat,IV)
    TS_u_mat <- cbind(TS_u_mat,u)


    finish <- finish +1
    setTxtProgressBar(pb, finish)

  }

  cat("\n")


  ###############################################
  #remove NA instruments
  sel_j <- NULL
  for(j in 1:n_IV){
    if(sum(is.na(TS_IV_mat[,j]))==0)sel_j<-c(sel_j,j)
  }

  #if there is NA elements in u_mat we delete it
  if(length(sel_j)>0){
    TS_IV_mat <- TS_IV_mat[,sel_j]
    TS_u_mat <- TS_u_mat[,sel_j]
  }

  ##################################################
  #now remove highly correlated IVs

  #IV.reduce <- IV_reduction(TS_IV_mat,crit_high_cor=0.9,maf_crit = (-1) )

  final_IV_mat <- TS_IV_mat
  final_u_mat  <- TS_u_mat

  #final_IV_mat <- IV.reduce$sel_snp
  #final_u_mat  <- TS_u_mat[,IV.reduce$id_snp]


  list(IV_mat= final_IV_mat, u_mat = final_u_mat, opt_lambda = opt_lambda,
       Pred_error_list =  Pred_error_list, G_pred_error_list = G_pred_error_list )

}

#' @title select IVs from a smooth_IV object (experimental function).
#' @description this function removes IVs with extreme low correlation and extreme high prediction error.
#' This is an experimental function to check how many redundant solutions are found in smooth.opt object.
#' @param smooth_IV: an object from smooth_CIV() function.
#' @param MR.data: data frame containing G,X,Z,Y.
#' @param ......: default values for other tuning parameters.
#' @return IV_mat: the final matrix of CIV instruments.
#' @return u_mat: the final CIV solutions of u. Each column is a distinct solution.
#' @examples
#' data(simulation)
#' G <- simulation$G
#' X <- simulation$X
#' Z <- simulation$Z
#' Y <- simulation$Y
#' smooth.opt <- smooth_CIV( G,X,Z,Y, k_folds = 10)
#' smooth.clean <- rm_outlier_IV(smooth.opt, simulation)
#' dim(smooth.clean$u_mat) #check how many solutions are different. It is probability much less than 100.
#' @export


rm_outlier_IV <- function(smooth_IV, MR.data, crit=0.9,sigma_min=0.01){

  X <- MR.data$X
  G <- MR.data$G
  Y <- MR.data$Y
  p<-ncol(G)
  IV_mat <- smooth_IV$IV_mat
  u_mat <- smooth_IV$u_mat

  cor_pen_list <- NULL
  er_list <- NULL
  G_er_list <- NULL

  for(j in 1:ncol(IV_mat)){

    temp_u <- smooth_IV$u_mat[,j]
    temp_IV <- IV_mat[,j]
    cor_pen_list <- c(cor_pen_list, cor(temp_IV,X) - p*smooth_IV$opt_lambda +
                        smooth_IV$opt_lambda *sum(exp(-temp_u^2/2/sigma_min^2))
    )


    MR.data <- MR.data
    MR.data$G <- temp_IV
    temp_tsls_fit <- TSLS_IV(MR.data)
    temp_Pred_error <- Y - cbind(1,X) %*% (temp_tsls_fit$coef)

    PG <- tcrossprod(temp_IV,temp_IV)/sum((temp_IV)^2)
    temp_G_pred_error <- PG %*% temp_Pred_error
    er_list <- c(er_list,sum(temp_Pred_error^2))
    G_er_list <- c(G_er_list, sum(temp_G_pred_error^2))
  }

  #we need to remove some columns of G according to a criteria(er_list or cor_list)

  #we identify the outliers in er_list and remove all these IVs
  id.er.outliers <- which(er_list > (  quantile(er_list,0.75) + 1.5*IQR(er_list) ) )

  #remove these IVs with too low correlation(penalized) with X
  id.cor.outliers <- which(cor_pen_list < (  quantile(cor_pen_list,0.25) - 1.5*IQR(cor_pen_list) ) )

  id_iv <- union(id.er.outliers,id.cor.outliers)

  if( length(id_iv)>0 ){
    IV_mat <- IV_mat[,-id_iv]
    u_mat <- u_mat[,-id_iv]
  }
  #######################################
  #now after obtain multiple IV we eliminate the replicated ones
  if(TRUE){
    select_u_mat  <- u_mat[,1]
    select_iv_mat <-IV_mat[,1]
    select_id <- id_iv[1]
    for(i in 2:ncol(IV_mat)){
      if( max(abs(cor(select_iv_mat,IV_mat[,i])))<crit){
        #print( c( dim(select_iv_mat),i) )
        select_iv_mat <- cbind(select_iv_mat,IV_mat[,i] )
        select_u_mat <- cbind(select_u_mat, u_mat[,i])
        select_id <- c(select_id,id_iv[i])
      }
    }
   u_mat <- select_u_mat
   IV_mat <- select_iv_mat
   id_iv <- select_id
  }

  list(IV_mat = IV_mat, u_mat = u_mat, id_iv=id_iv, er_list =er_list,
       G_er_list = G_er_list, cor_pen_list = cor_pen_list )



}

#' @title Two stage least square method.
#' @description This function implement ordinary two stage least square regression and provide variance estimation (if requested).
#' @param MR.data: data frame containing G,X,Z,Y.
#' @param Fstats: return F-statistics or not. If multiple phenotypes (X) are used, Pillai statistics will be used instead.
#' @param var_cal: return variance estimation or not.
#' @return coef: the causal effect estimation \eqn{\beta}.
#' @return var: the variance estimation of \eqn{\beta}. if var_cal=TRUE.
#' @return stats: F-statistics (or Pillai statistics). if Fstats=TRUE.
#' @return pvalue: the pvalue of F-statistics. if Fstats=TRUE.
#' @examples
#' data(simulation)
#' TSLS_IV(simulation,Fstats=TRUE,var_cal=TRUE)
#' @export


TSLS_IV <- function(MR.data,Fstats=FALSE,var_cal=FALSE){


  MR.data$Y -> Y

  MR.data$X -> X

  MR.data$G -> G

  if(length(ncol(MR.data$X))==0)X <- matrix(MR.data$X,ncol=1)


  p <- ncol(G)

  if(sum(p)==0){p <- 1}

  n <- length(Y)


  #the method in sem package
  #tsls_fit <- tsls(y~X, ~G, data=MR.data)

  #the naive implementation
  x_fit<- lm(X~G,data=MR.data)$fitted.values

  x_fit<- lm((MR.data$X)~(MR.data$G) )$fitted.values


  y_fit<- lm(MR.data$Y ~ x_fit)

  #calculate the beta_iv
  beta_tsls <- y_fit$coefficients

  var_tsls <- 0

  if(var_cal==TRUE){
    #calcuate the asymtotic variance of this estimator
    Pz <- G%*% solve(t(G)%*%G) %*% t(G)

    xzzz_inv <- t(X)%*%G %*% solve(t(G)%*%G)
    xPzx <- solve(t(X)%*%Pz%*%X)

    tsls_fit <- tsls(Y~X, ~G, data=MR.data)
    resid <- tsls_fit$residuals


    current_mat <- matrix(0,p,p)

    if(p>1){
      for(i in 1:n){
        current_mat <- current_mat + ((resid[i])^2)*G[i,]%*%t(G[i,])
      }
    }

    else{
      for(i in 1:n){
        current_mat <- current_mat + ((resid[i])^2)*(G[i])^2
      }


    }

    S <- current_mat / n

    var_tsls <- n*xPzx%*%( xzzz_inv %*% S %*% t(xzzz_inv)) %*% xPzx
  }

  if(Fstats==TRUE){
    #now give out the 1st stage F statistics (strong instrument or not) and its pvalue

    if(ncol(X)==1){
      lmfit <- lm(X~G)
      stats <- summary(lmfit)$fstatistic[1]
      pvalue <- (anova(lmfit)$`Pr(>F)`)[1]
    }

    else{

      MVR <- manova(X~G)

      #record Pillai statistics and p value, and association with y
      stats <- summary(MVR)$stats[1,"Pillai"]
      pvalue <- summary(MVR)$stats[1,"Pr(>F)"]

    }

    #coef is the coefficient estimator of y|X. var is the variance of this estimator.
    #stat and pvalue is the stage I assiciation test result.
    list(coef=beta_tsls,var=var_tsls,stats = stats, pvalue=pvalue)
  }

  else{
    list(coef=beta_tsls,var=var_tsls)
  }
}

#' @title cross-validated Allele score method.
#' @description This function implement Allele score methods with cross-validation
#' in the way Stephen Burgess suggested in the Allele score methods paper.
#' @param MR.data: data frame containing G,X,Z,Y.
#' @param n_folds: the number of folds for cross-validation.
#' @return weights: the weights for allele score across folds. Each row is a weight vector
#' corresponding to a specific fold.
#' @return allele_score: The cross-validated Allele score, which would be used as the new instruments
#' in MR analysis.
#' @return beta_est: the causal effect estimation of \eqn{\beta}.
#' @examples
#' data(simulation)
#' allele.score <- allele(simulation,n_folds=10)
#' @export


allele <- function(MR.data, n_folds = 10 ){

  G <- MR.data$G
  X <- MR.data$X
  Z <- MR.data$Z
  Y <- MR.data$Y

  if(length(ncol(MR.data$X))==0)X <- matrix(MR.data$X,ncol=1)
  n <- length(Y)
  p <- ncol(G)

  #allele_score is used to save the calculated allele scores
  allele_score <- X
  #weights is a matrix where each row contains weight from each folds
  weights <- NULL
  flds <- createFolds(c(1:n), k = n_folds, list = TRUE, returnTrain = FALSE)

  for(cv_folds in 1:n_folds){

    apply_index <- flds[[cv_folds]]
    train_index <- c(1:n)[-apply_index]

    if(ncol(X)==1){X_train <- as.matrix(X[train_index,1],ncol=1)
    X_apply <- as.matrix(X[apply_index,1],ncol=1)
    }
    else{ X_train <- X[train_index,]   }

    G_train <- G[train_index,]
    G_apply  <- G[apply_index,]

    coef_train <- lm(X_train~G_train)$coefficients

    coef_train[which(is.na(coef_train))] <- 0

    if(ncol(X)==1){
      if(ncol(G)>1){
        allele_score_apply <- G_apply%*%coef_train[-1]}
      else{
        allele_score_apply <- G_apply * coef_train[-1]
      }
      allele_score[apply_index] <- allele_score_apply
    }


    if(ncol(X)>1){
      allele_score_apply <- G_apply%*%coef_train[-1,]
      allele_score[apply_index,] <- allele_score_apply
    }


    weights <- rbind(weights,coef_train)
  }

  #now we use the allele score to do causal effect estimation
  temp_TSLS_data <- MR.data
  temp_TSLS_data$G <- allele_score

  est_TSLS <- TSLS_IV(temp_TSLS_data)
  beta_est <- est_TSLS$coef[-1]


  list(beta_est = beta_est, weights = weights, allele_score = allele_score)



}

#' @title simple linear regression pvalues (internal function.)
#' @description univaraite t-test pvalues for a regression.
#' @param modelobject: a regression object.
#' @return p: pvalue.
#' @export
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#' @title univariate T-test p-values.
#' @description Given response $Y$ and a set of features $X$, this
#' function obtains univariate T-test p-values for each of the feature for selection purpose.
#' @param Y: response variable. \eqn{n \times 1}.
#' @param X: the independent features.\eqn{n \times p}.
#' @return pvalue_list: the list of Pvalues for feature selection based on univariate T-test.
#' @examples
#' data(simulation)
#' p.values <- lmPvalue(simulation$X, simulation$G)
#' @export
lmPvalue <- function(Y,X){

  pvalue_list <- NULL

  for(j in 1:ncol(X)){
    fit <- (lm(Y~X[,j]))
    pvalue_list <- c( pvalue_list, lmp(fit) )
  }
  pvalue_list
}

