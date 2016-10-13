#' Regularized graphical model estimation
#'
#' \code{greg} calculate the regularized graphical model estimation using lasso, scad and adaptive lasso penalties. It report the results in the form of roc results for each method.
#'
#' @export
#' @importFrom parallel mcmapply
#' @param z n * p dimensional matrix
#' @param A p * p true graph
#' @param eps a tolerence level for thresholding
#' @param rholist a sequence of penalty parameters
#' @param gamma the adaptive lasso penalty parameter
#' @param trace whether to trace to estimation process.
#' @return a list.
#' \item{roc.lasso}{roc results for lasso}
#' \item{roc.scad}{roc results for scad}
#' \item{roc.alasso}{roc results for adaptive lasso}
#' @seealso \code{\link{pgraph}}, \code{\link{roc}}, \code{\link{projcov}}
#' @examples
#' set.seed(0)
#' p = 20;
#' n = 300;
#' tmp=runif(p-1,1,3)
#' s=c(0,cumsum(tmp));
#' s1=matrix(s,p,p)
#' cov.mat.true=exp(-abs(s1-t(s1)))
#' prec.mat.true=solve(cov.mat.true);
#' a=matrix(rnorm(p*n),n,p)
#' data.sa=a%*%chol(cov.mat.true);
#' true.graph = outer(1:p,1:p,f<-function(x,y){(abs(x-y)==1)})
#' greg.fit = greg(data.sa, true.graph)
#' auc.lasso = sum(diff(greg.fit$roc.lasso[,1])*greg.fit$roc.lasso[-1,2])
#' auc.alasso = sum(diff(greg.fit$roc.alasso[,1])*greg.fit$roc.alasso[-1,2])
#' auc.scad = sum(diff(greg.fit$roc.scad[,1])*greg.fit$roc.scad[-1,2])
#' auc.lasso
#' auc.alasso
#' auc.scad
greg <- function(z, A, eps = 1e-15, rholist = NULL, gamma = 0.5, trace = FALSE){
 diag(A) = 0 ###true graph has diagonal 0
# n = nrow(z)
 p = ncol(z)
 s = cov(z)
 if (is.null(rholist)) {
   rholist = exp(seq(log(max(abs(s))/1000), log(max(abs(s))), length = 50))
 }
 rholist = sort(rholist, decreasing = T)
 n.rho = length(rholist)
 FP.lasso = FN.lasso = FP.alasso = FN.alasso = FP.scad = FN.scad = rep(0,n.rho)
 for (i in 1:n.rho) {
   rho = rholist[i]
   fit.lasso = glasso::glasso(s, rho)
   wi.lasso = fit.lasso$wi
   Ahat.lasso = (abs(wi.lasso) > eps)*1
   diag(Ahat.lasso) = 0
   FP.lasso[i] = sum(Ahat.lasso > A)
   FN.lasso[i] = sum(Ahat.lasso < A)

   rhomat = rho/p/2 * matrix(1,p,p)/(pmax(abs(wi.lasso) ^ gamma,1e-5));
   fit.alasso = glasso::glasso(s, rhomat);
   wi.alasso = fit.alasso$wi
   Ahat.alasso = (abs(wi.alasso) > eps)*1
   diag(Ahat.alasso) = 0
   FP.alasso[i] = sum(Ahat.alasso > A)
   FN.alasso[i] = sum(Ahat.alasso < A)

   wi.scad = wi.lasso;
   epsi = 1
   count = 1
   while (epsi > 1e-4) {
     if (epsi < 1e-3 & count > 20) break
     count = count + 1;
     wi.scad.old = wi.scad;
     rhomat = scadrightderv(abs(wi.scad.old),3.7,rho);
     fit.scad = glasso::glasso(s, rhomat)
     wi.scad = fit.scad$wi;
     epsi = mean(abs(wi.scad - wi.scad.old));
     if ( count > 50)
     {
       warning("scad iteration does not converge")
       break
     }
   }

   Ahat.scad = (abs(wi.scad) > eps)*1
   diag(Ahat.scad) = 0
   FP.scad[i] = sum(Ahat.scad > A)
   FN.scad[i] = sum(Ahat.scad < A)

   if (trace) cat('i=', i, ' ', FP.lasso[i], FN.lasso[i], FP.alasso[i], FN.alasso[i], FP.scad[i], FN.scad[i],'\n')
 }

 true1 = which(A > 0)
 true0 = which(A == 0)
 n1 = length(true1)
 n0 = length(true0)

 roc.lasso = cbind(fpr = FP.lasso/(n0 - p), tpr = 1 - FN.lasso/n1)

 roc.scad = cbind(fpr = FP.scad/(n0 - p), tpr = 1 - FN.scad/n1)

 roc.alasso = cbind(fpr = FP.alasso/(n0 - p), tpr = 1 - FN.alasso/n1)


return(list(roc.lasso = roc.lasso,roc.scad = roc.scad, roc.alasso = roc.alasso))

}
