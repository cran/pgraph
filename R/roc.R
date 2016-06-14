#' Compute the Projected Graph
#'
#' \code{roc} calculate the fpr and tpr for the roc curve
#'
#' @export
#' @param a p * p estimated graph
#' @param a0 p * p true graph
#' @return a list.
#' \item{tpr}{tpr sequence}
#' \item{fpr}{fpr sequence}
#' @examples
#' require(splines)
#' p = 20;
#' n = 50;
#' tmp=runif(p-1,1,3)
#' s=c(0,cumsum(tmp));
#' s1=matrix(s,p,p)
#' cov.mat.true=exp(-abs(s1-t(s1)))
#' prec.mat.true=solve(cov.mat.true);
#' a=matrix(rnorm(p*n),n,p)
#' data.sa=a%*%chol(cov.mat.true);
#' true.graph = outer(1:p,1:p,f<-function(x,y){(abs(x-y)==1)})
#' #fit = greg(data.sa, true.graph)
#' #plot(fit$roc.lasso[,1],fit$roc.lasso[,2],type='l')
#' #lines(fit$roc.alasso[,1],fit$roc.alasso[,2],col=2)
#' #lines(fit$roc.scad[,1],fit$roc.scad[,2],col=3)
#' #methodlist = c("lasso","sam")
#' #fit = vector(mode="list", length=2)
#' #info = vector(mode="list", length=2)
#' #auc = NULL
#' #plot.new()
#' #for(i in 1:2){
#' #method = methodlist[i]
#' #fit[[i]] = pgraph(data.sa, method = method)
#' #info[[i]] = roc(fit[[i]]$statmat.pearson, true.graph)
#' #auc[i] = sum(-diff(info[[i]][,1])*info[[i]][-1,2])
#' #lines(info[[i]], xlab='FPR',ylab='TPR', type='s', col=i+3)
#' #cat(method, ': auc=', auc[i],'\n')
#' #}
#' #auc
roc <- function(a, a0){
 cutoff = unique(sort(as.vector(a)))
 n.cutoff = length(cutoff)
 fpr = tpr = rep(0,n.cutoff)
 true1 = which(a0>0)
 true0 = which(a0==0)
 n = nrow(a)
 for(i in 1:n.cutoff){
   est1 = which(a>cutoff[i])
   est0 = which(a<=cutoff[i])
   tpr[i] = length(intersect(est1,true1))/length(true1)
   fpr[i] = length(setdiff(est1,true1))/(length(true0)-n)##ignore the digonal elements
 }
 return(cbind(fpr=fpr,tpr=tpr))


}
