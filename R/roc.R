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
