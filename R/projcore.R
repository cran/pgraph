#' Calculate the Projected matrix given factors
#'
#' \code{projcore} calculate the projected matrix given
#' factors.
#'
#' @export
#' @param x first vector
#' @param b factor matrix
#' @param method projection method. Default = 'linear'.
#' @param one.SE whether to use the 1se rule for glmnet. Default = TRUE.
#' @param refit whether to refit the selected model. Default = TRUE.
#' @param randSeed the random seed for the program. Default = 0.
#' @return
#' \item{eps}{the residual matrix after projection}
#' @seealso \code{\link{greg}}, \code{\link{roc}}, \code{\link{pgraph}}
projcore <- function(x, b, method = c("lasso","sam"), one.SE = TRUE, refit = TRUE, randSeed = 0){
  method = match.arg(method)
  # cor = match.arg(cor)
  set.seed(randSeed)
  x = as.matrix(x)
  p = ncol(x)
  resi = x
  for (j in 1:p) {
    if(method == 'lasso'){
      xfit = glmnet::cv.glmnet(b,x[,j])
      if(one.SE){
        if(refit){
          bi = cbind(1,b)
          index = which(as.vector(coef(xfit, s = xfit$lambda.1se)!=0))
          fit = lm.fit(as.matrix(bi[,index]), x[,j])
          resi[,j] = fit$residuals
        } else{
          resi[,j] = x[,j] - predict(xfit, b, xfit$lambda.1se)
        }
      } else{
        if(refit){
          bi = cbind(1,b)
          index = which(as.vector(coef(xfit, s = xfit$lambda.min)!=0))
          fit = lm.fit(as.matrix(bi[,index]), x[,j])
          resi[,j] = fit$residuals
        } else{
          resi[,j] = x[,j] - predict(xfit, b, xfit$lambda.min)
        }
      }
    } else if(method == 'sam'){
      xfit = cv.samQL(b,x[,j])
      resi[,j] = x[,j] - predict(xfit$sam.fit, b)$values[,xfit$lambda.min]
    }
  }


  return(resi = resi)
}
