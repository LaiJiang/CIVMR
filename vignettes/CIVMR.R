## ----  results='asis',message=FALSE, warning=FALSE-----------------------
      library(CIVMR)
      data(simulation)
      G <- simulation$G
      colnames(G) <- paste("rs", seq(1:9), sep="")
      rownames(G) <- paste("sub", seq(1:nrow(G)), sep="")
      knitr::kable(G[1:5,])

## ----out.width='\\textwidth', fig.width = 5, fig.height = 5, fig.align='center',message=FALSE, warning=FALSE----
library(CIVMR)
data(list=ADNI)
list_mat <- NULL
list_mat <- rbind(list_mat, lmPvalue(ADNI$X,ADNI$G) )
list_mat <- rbind(list_mat, lmPvalue(ADNI$Z[,1],ADNI$G) )
list_mat <- rbind(list_mat, lmPvalue(ADNI$Z[,2],ADNI$G) )
list_mat <- rbind(list_mat, lmPvalue(ADNI$Z[,3],ADNI$G) )
list_mat <- log(list_mat)

plot(list_mat[1,],col="red", ylab ="log-Pvalue",  pch =20, ylim=c(0,min(list_mat)),
     xlab="20 SNPs associated with Amyloid beta 1-42 in previous studies")
text(5, -20, "Amyloid beta" ,col="red" )
text(5, -30, "T-Tau" ,col="blue" )
text(5, -40, "P-Tau" ,col="black" )
text(5, -50, "FDG: SUVR_global" ,col="green" )
points(list_mat[2,], col="blue", pch =20 )
points(list_mat[3,], col="black", pch =20)
points(list_mat[4,], col="green", pch =20)
abline(h=log(0.05))
abline(h=log(5*10^(-8)))


## ----  results='asis',message=FALSE, warning=FALSE,results='hide'--------
library(AER)
fm <- ivreg(Y ~ X | G, data = ADNI)
summary(fm, vcov = sandwich, df = Inf, diagnostics = TRUE)

##Call:
##ivreg(formula = Y ~ X | G, data = ADNI)

##Residuals:
##    Min      1Q  Median      3Q     Max 
##-0.9314 -0.4743  0.1744  0.3129  0.7423 

##Coefficients:
##              Estimate Std. Error z value Pr(>|z|)    
##(Intercept) -1.513e-17  1.958e-02    0.00        1    
##X           -2.501e-03  5.777e-04   -4.33 1.49e-05 ***

##Diagnostic tests:
##                 df1 df2 statistic p-value    
##Weak instruments  20 470    17.375  <2e-16 ***
##Wu-Hausman         1 488     0.001   0.976    
##Sargan            19  NA    23.312   0.224    
##---
##Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Residual standard error: 0.4347 on Inf degrees of freedom
##Multiple R-Squared: 0.09217,	Adjusted R-squared: 0.09031 
##Wald test: 18.75 on 1 DF,  p-value: 1.494e-05 

## ----  results='asis',message=FALSE, warning=FALSE,results='hide'--------

G <- lm(G~Z,data=ADNI)$residuals
X <- lm(X~Z,data=ADNI)$residuals
Y <- lm(Y~Z,data=ADNI)$residuals
ADNI.adj <- list(G=G,X=X,Y=Y)

ADNI.mul <- list(G=ADNI$G, X=cbind(ADNI$X,ADNI$Z),Y=ADNI$Y)

#2SLS method (naive)
tsls.naive <- TSLS_IV(ADNI, Fstats=TRUE, var_cal=TRUE)
#2SLS method (adjusted)
tsls.adj <- TSLS_IV(ADNI.adj, Fstats=TRUE, var_cal=TRUE)
#2SLS method (multiple regression)
ADNI.mul <- list(G=ADNI$G, X=cbind(ADNI$X,ADNI$Z),Y=ADNI$Y)
tsls.mul <- TSLS_IV(ADNI.mul, Fstats=TRUE, var_cal=TRUE)

## ----  results='asis',message=FALSE, warning=FALSE,results='hide'--------

#1. naive allele score method.
allele.est <- allele(ADNI)

#2. allele score method while adjusting for pleiotropic phenotypes Z
allele.adj <- allele(ADNI.adj)

## ----  results='asis',message=FALSE, warning=FALSE,results='hide'--------
#this step is not mandatory
sel.G <- IV_reduction(ADNI$G, crit_high_cor = 0.9)$sel_snp
#then calculate the CIV instrument
civ.G <- CIV(list(G=sel.G,X=ADNI$X,Z=ADNI$Z,Y=ADNI$Y) )$CIV
#then run the TSLS estimation method
civ.est <- TSLS_IV(list(G=civ.G, X = ADNI$X, Y=ADNI$Y, y=ADNI$Y), Fstats=TRUE, var_cal=TRUE)

#2. then the bootstrapped CIV instrument
civ.boot.G <- sel.G %*% (boot_CIV(list(G=sel.G,X=ADNI$X,Z=ADNI$Z,Y=ADNI$Y) )$boots.cor.u)
#civ.boot estimation.
civ.boot.est <- TSLS_IV(list(G= civ.boot.G, X = ADNI$X, Y=ADNI$Y, y=ADNI$Y), Fstats=TRUE, var_cal=TRUE)

#3. then the croos-validated CIV instrument.
civ.cv.G <- cv_CIV(list(G=sel.G,X=ADNI$X,Z=ADNI$Z,Y=ADNI$Y))

## ----  results='asis',message=FALSE, warning=FALSE,results='hide'--------
#Find CIV_smooth solution.
civ.smooth <- smooth_CIV(ADNI$G,ADNI$X,ADNI$Z,ADNI$Y)

#this step is not mandatory
civ.iv <- IV_reduction(civ.smooth$IV_mat, crit_high_cor = 0.9)

#now we obtain the causal effect estimation.
civ.smooth.est <- TSLS_IV(list(G=civ.iv$sel_snp, X = ADNI$X, Y=ADNI$Y, y=ADNI$Y), Fstats=TRUE, var_cal=TRUE)


## ---- fig.show='hold'----------------------------------------------------
plot(civ.smooth$u_mat[,civ.iv$id_snp[1]],xlab="Solution 1",ylab="weights")
plot(civ.smooth$u_mat[,civ.iv$id_snp[2]],xlab="Solution 2",ylab="weights")

## ---- echo=FALSE, results='asis'-----------------------------------------

ce.table <- c(tsls.naive$coef[-1], tsls.adj$coef[-1], tsls.mul$coef[2], 
              allele.est$beta_est, allele.adj$beta_est, 
              civ.est$coef[-1], civ.boot.est$coef[-1], civ.cv.G$beta_est,
              civ.smooth.est$coef[-1]
              )

var.table <- c(tsls.naive$var[1], tsls.adj$var[1], tsls.mul$var[1],
               NA,NA,
               civ.est$var[1], civ.boot.est$var[1], NA,
               civ.smooth.est$var[1])

result.table <- rbind(ce.table, var.table)

colnames(result.table) <- c("2SLS.naive","2SLS.adj","2SLS.mul",
                            "Allele.naive","Allele.adj",
                            "CIV","CIV_boot","CIV_cv",
                            "CIV_smooth")
rownames(result.table) <- c("causal effect estimation","variance estimation")
knitr::kable(head(result.table, 9))

