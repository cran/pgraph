#' Calculate the Conditional Dependency Graph
#'
#' \code{pgraph} calculate the conditional dependency graph (with/without external factors) via projection using lasso or sparse additive model.
#'
#' @export
# @importFrom parallel mcmapply
#' @param z n * p dimensional matrix
#' @param f n * q factor matrix. Default = 'NULL'.
#' @param method projection method. Default = 'linear'.
#' @param cond whether to create a conditional graph or unconditional graph.
#' Default = TRUE. If cond = FALSE, f must be provided.
#' @param randSeed the random seed for the program. Default = 0.
#' @param trace whether to trace to estimation process.
#' @return a list to be used to calculate the ROC curve.
#' \item{statmat.pearson}{matrix with pearson correlation test}
#' \item{statmat.dcov}{matrix with distance covariance test}
#' @seealso \code{\link{greg}}, \code{\link{roc}}, \code{\link{projcov}}
#' @examples
#' library(splines)
#' set.seed(0)
#' p = 5;
#' n = 100;
#' tmp=runif(p-1,1,3)
#' s=c(0,cumsum(tmp));
#' s1=matrix(s,p,p)
#' cov.mat.true=exp(-abs(s1-t(s1)))
#' prec.mat.true=solve(cov.mat.true);
#' a=matrix(rnorm(p*n),n,p)
#' data.sa=a%*%chol(cov.mat.true);
#' true.graph = outer(1:p,1:p,f<-function(x,y){(abs(x-y)==1)})
#' methodlist = c("lasso","sam")
#' fit = vector(mode="list", length=2)
#' info = vector(mode="list", length=2)
#' auc = NULL
#' for(i in 1:2){
#' method = methodlist[i]
#' fit[[i]] = pgraph(data.sa, method = method)
#' info[[i]] = roc(fit[[i]]$statmat.pearson, true.graph)
#' auc[i] = sum(-diff(info[[i]][,1])*info[[i]][-1,2])
#'   cat(method, ': auc=', auc[i],'\n')
#' }
pgraph <- function(z, f = NULL, method = c("lasso","sam"), cond = TRUE, randSeed = 0, trace = FALSE){
 method = match.arg(method)
 set.seed(randSeed)
 p = ncol(z)
 statmat.pearson = matrix(0,p,p)
 statmat.dcov = matrix(0,p,p)
 for(i in 1:(p-1))
  for(j in (i+1):p){
    if(trace){
      cat('i=',i,'j=',j,'\n')
    }

x = z[,i]
y = z[,j]
if(cond) {
  b = cbind(z[,-c(i,j)],f)
} else {
  if(is.null(f)) {
    stop('f must be provided when cond is FALSE.')
  }
  b = f
}
fit = projcov(x,y,b,method, randSeed = randSeed)
statmat.pearson[i,j] = fit$test.pearson
statmat.dcov[i,j] = fit$test.dcov
  }
 statmat.pearson = statmat.pearson+t(statmat.pearson)
 statmat.dcov = statmat.dcov+t(statmat.dcov)
return(list(statmat.pearson = statmat.pearson,statmat.dcov = statmat.dcov))

}
