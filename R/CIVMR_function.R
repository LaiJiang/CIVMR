
#' @title linear algebra decompositions for CIV
#' @description This function implements linear algebra steps to acquire necessary matrices for CIV construction.
#' @param G: original instruments with dimension nXp.
#' @param X: phenotype of interest. dimension nXk.
#' @param Z: possible pleiotropic phenotypes which have been measured. dimension nXr.
#' @keywords Cholesky decomposition
#' @return A list of matrices which will be called by solv_pcc() and pcc_IV().
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





#' @title Find a unique solution of CIV.
#' @description This function find a unique solutin to the constrained instrument problem given matrix A and B.
#' @param A: matrix given by LA_decomposition().
#' @param B: matrix given by LA_decomposition().
#' @keywords Cholesky decomposition
#' @return c: solution to the constrained maximization problem.
#' @return max_value: the maximized correlation value.
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



#' @title multiple orthogonal CIV solutions (if applicable)
#' @description This function find multiple CIV solutions that are orthogonal to each other if called. Only the first one achive the maximum correlation.
#' @param A: matrix given by LA_decomposition().
#' @param B: matrix given by LA_decomposition().
#' @param G: original instruments.
#' @return u_max: the solution of u that would maximize the constrained correlation problem.
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


#' @title cross-validated CIV.
#' @description This function generate a CIV using cross-validation.
#'              Specifically, for a predefined fold the CIV is calculated using
#'               all samples except this fold, then the CIV solution is applied to the samples in this
#'               fold to obtain CIV instrumentable variable. In this way the correlation between samples
#'               are expected to be reduced.
#' @param temp_TSLS_data: a data frame containing G,X,Z,Y.
#' @param n_folds: number of folds for cross-validation.
#' @return civ.IV: cross-validated CIV instrument.
#' @return beta_est: causal effect estimation of X on Y using civ.IV.
#' @export
cv_CIV <- function(temp_TSLS_data, n_folds = 10 ){

  G <- temp_TSLS_data$G
  X <- temp_TSLS_data$X
  Z <- temp_TSLS_data$Z
  Y <- temp_TSLS_data$Y

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
  new_TSLS_data <- temp_TSLS_data
  new_TSLS_data$G <- civ_score

  est_TSLS <- TSLS_IV(new_TSLS_data)
  beta_est <- est_TSLS$coef[-1]


  list(beta_est = beta_est, weights = weights, civ.IV = civ_score)
}


#' @title bootstrapped CIV.
#' @description This function generate a bootstrapped CIV w/wo correction.
#'              Specifically, for a bootstrap sample we can generate civ solution u. The boostrap corrected
#'              solution u is obtained as the global solution u - ( bootrapped average u - global u).
#' @param temp_TSLS_data: a data frame containing G,X,Z,Y.
#' @param n_boots: number of bootstrap samples.
#' @return boots.u: bootstrapped CIV solution u (without correction).
#' @return boots.cor.u: bootstrap corrected solution of u.
#' @export


boot_CIV <- function(temp_TSLS_data, n_boots = 10 ){

  X <- temp_TSLS_data$X
  G <- temp_TSLS_data$G
  Z <- temp_TSLS_data$Z
  Y <- temp_TSLS_data$Y

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
    cor_z <- cor(boot_Z)
    diag(cor_g) <- diag(cor_z) <- 0

    #we keep re-sampling if the previous sampled Gs are highly correlated
    while( (is.na(max(cor_g))==1)||(max(abs(cor_g))>0.999)||(max(abs(cor_z))>0.999) ){

      obs_id <- sample(1:n, size=n, replace=TRUE)

      boot_X <- X[obs_id,]
      boot_Z <- Z[obs_id,]
      boot_G <- G[obs_id,]
      boot_Y <- Y[obs_id]

      cor_g <- cor(boot_G)
      cor_z <- cor(boot_Z)
      diag(cor_g) <- diag(cor_z) <- 0
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


  deco <- LA_decomposition(temp_TSLS_data$G,temp_TSLS_data$X,temp_TSLS_data$Z)
  pcc_iv <- pcc_IV(deco$A,deco$B,temp_TSLS_data$G,
                   deco$inv_GG_square,no_IV = 1)
  #weight
  civ.global.u <- pcc_iv$u_max

  boots.u <- apply(civ.boot.mat,1,mean)

  boots.cor.u <- 2*civ.global.u - boots.u

  list(boots.u = boots.u, boots.cor.u = boots.cor.u)
}



#' @title SNP pre-processing.
#' @description This function remove highly correlated SNPs. It also
#' calculate MAF for each snp and delete rare snps with 0.01 MAF

#' @param snp_matrix: SNP matrix with dimension n * p.
#' @param crit_high_cor: criteria to choose highly correlated SNPs.
#' @param maf_crit: criteria to choose rare SNPs.
#' @return sel_snp: the selected SNPs.
#' @return id_snp: the ids (columns) of selected SNPs in the original SNP matrix.
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
#' @description This function remove highly correlated IVs.

#' @param snp_matrix: IV matrix with dimension n * p.
#' @param crit_high_cor: criteria to choose highly correlated SNPs.
#' @return sel_snp: the selected IVs.
#' @return id_snp: the ids (columns) of selected IVs in the original IV matrix.
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

#' @title CIV_smooth with a specified lambda.
#' @description This function finds a CIV_smooth solution of u given a lambda value.
#' @param initial: the initial value for updating u.
#' @param G: SNP matrix of n*p.
#' @param X: phenotype of interest.
#' @param lambda: a given value (must be specified) for regularization parameter.
#' @param ......: default values for other tuning parameters.
#' @return mat_u: the trace of all updated iterations of u.
#' @return opt_solution: the converged final solution of u.
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
        pcc_iv <- pcc_IV(A,B,G, inv_GG_square)
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

#' @title CIV_smooth solution with cross-validation.
#' @description This function first find the optimal lambda value according
#' to projected prediction error with cross-validation. Then for a given lambda value multiple intial
#' points are used to explore potentially multiple modes.
#' @param initial: the initial value for updating u.
#' @param G: SNP matrix of n*p.
#' @param X: phenotype of interest.
#' @param lambda_list: a list of values for regularization parameter lambda.
#' @param k_folds: number of folds for cross-validation.
#' @param n_IV: the number of initial points chosen to explore potential multiple modes. The converged
#' solutions will be screened to delete redundant solutions. So the final solutions will be less than n_IV.
#' @param ......: default values for other tuning parameters.
#' @return IV_mat: the final matrix of CIV instruments.
#' @return u_mat: the final CIV solutions of u.
#' @return G_pred_error_list: the projected prediction error.
#' @export
smooth_opt_IV <- function(lambda_list = NULL, k_folds =2,
                          sigma_min = 0.01, sigma_up =0.5,
                          stepsize=0.1, conv_iters = 5,
                          stepsize_last = 0.0001,
                          last_conv_iters = 2000,
                          null_space, G,X,Z,Y, method_lambda ="er",
                          n_IV = 2){

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
      temp_TSLS_data <- data.frame(Y_test)
      temp_TSLS_data$G <- IV_test
      temp_TSLS_data$X <- X_test
      temp_TSLS_data$y <- Y_test
      cv_tsls_fit <- TSLS_IV(temp_TSLS_data)
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


#' @title select IVs according to their training prediction error (if necessary).
#' @description This function selects IVs according to their training prediction error.(Y-X *beta)
#' @param smooth_IV: object from smooth_opt_IV() function.
#' @param TSLS_data: data frame containing G,X,Z,Y.
#' @param ......: default values for other tuning parameters.
#' @return IV_mat: the final matrix of CIV instruments.
#' @return u_mat: the final CIV solutions of u.
#' @export

select_IV <- function(smooth_IV, TSLS_data, crit=0.9,sigma_min=0.01){

  X <- TSLS_data$X
  G <- TSLS_data$G
  Y <- TSLS_data$Y
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


    temp_TSLS_data <- TSLS_data
    temp_TSLS_data$G <- temp_IV
    temp_tsls_fit <- TSLS_IV(temp_TSLS_data)
    temp_Pred_error <- Y - cbind(1,X) %*% (temp_tsls_fit$coef)

    PG <- tcrossprod(temp_IV,temp_IV)/sum((temp_IV)^2)
    temp_G_pred_error <- PG %*% temp_Pred_error
    er_list <- c(er_list,sum(temp_Pred_error^2))
    G_er_list <- c(G_er_list, sum(temp_G_pred_error^2))
  }

  #we need to choose the columns of G according to a criteria(er_list or cor_list)
  median_er <- median(er_list)
  min_er <- min(er_list)
  up_bound <-  median_er + 1*(median_er-min_er)
  below_bound <- median_er - 1*(median_er-min_er)
  #id_iv <- which( (er_list<up_bound)*(er_list > below_bound)==1    )
  id_er <- which(er_list < up_bound)

  median_cor <- median(cor_pen_list)
  max_cor <- max(cor_pen_list)
  #id_iv <- which( (er_list<up_bound)*(er_list > below_bound)==1    )
  id_cor <- which(cor_pen_list > quantile(cor_pen_list,0.5) )

  id_iv <- intersect(id_er,id_cor)
  #id_iv <- id_cor

  IV_mat <- IV_mat[,id_iv]
  u_mat <- u_mat[,id_iv]
  #######################################
  #now after obtain multiple IV we eliminate the replicated ones
  if(FALSE){
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

  }

  list(IV_mat = IV_mat, u_mat = u_mat, id_iv=id_iv, er_list =er_list,
       G_er_list = G_er_list, cor_pen_list = cor_pen_list )



}


#' @title select IVs according to their training prediction error and correlation with X (if necessary).
#' @description #this function removes IVs with extreme low correlation and extreme high prediction error.
#' @param smooth_IV: object from smooth_opt_IV() function.
#' @param TSLS_data: data frame containing G,X,Z,Y.
#' @param ......: default values for other tuning parameters.
#' @return IV_mat: the final matrix of CIV instruments.
#' @return u_mat: the final CIV solutions of u.
#' @export


rm_outlier_IV <- function(smooth_IV, TSLS_data, crit=0.9,sigma_min=0.01){

  X <- TSLS_data$X
  G <- TSLS_data$G
  Y <- TSLS_data$Y
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


    temp_TSLS_data <- TSLS_data
    temp_TSLS_data$G <- temp_IV
    temp_tsls_fit <- TSLS_IV(temp_TSLS_data)
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
  if(FALSE){
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

  }

  list(IV_mat = IV_mat, u_mat = u_mat, id_iv=id_iv, er_list =er_list,
       G_er_list = G_er_list, cor_pen_list = cor_pen_list )



}




