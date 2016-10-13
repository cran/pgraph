cv.samQL <- function(x, y, fold = 5){
  n = length(y)
  n.fold = ceiling(n/fold); # number of observations in each fold
  set.seed(0)
  ind = sample(1:n)
  error = rep(0,30)
  for(i in 1:fold){
    cind = (i - 1) * n.fold + 1:n.fold
    cind = intersect(1:n, cind)
    test.ind = ind[cind]
    train.ind = setdiff(ind,test.ind)
    fit = SAM::samQL(x[train.ind,],y[train.ind])
    ypred = predict(fit, x[test.ind,])
    error = error + apply(ypred$values,2,f<-function(u){sum((u-y[test.ind])^2)})
  }
  lambda.min = which.min(error)
  sam.fit = SAM::samQL(x,y)
  return(list=list(sam.fit=sam.fit,lambda.min=lambda.min,cv.error=error))
}

scadrightderv <- function(lamhat, a, lam)
{
  pmax(lam*((lamhat<=lam)+pmax(a*lam-lamhat,0)*(lamhat>lam)/(a-1)/lam),1e-10)
}

