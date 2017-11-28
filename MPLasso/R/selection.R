#####################################################################
# Graphical Lasso optimal lambda selection
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

selection <- function(sample, type = "MPLasso", prior = 1, priorType = TRUE, interactionFlag = FALSE, precision=FALSE, precisionRatio = 0.5, priorInfo){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);
  
  lambdaMin = 0.01;
  lambdaMax = 1;
  interval = 0.01;
  numLambda = (lambdaMax - lambdaMin)/interval + 1;
  BICResult = matrix(0, numLambda, 1);
  
  for (i in 1:numLambda){
    result <- MPLasso(sample, lambdaMin = (lambdaMin+(i-1)*interval), lambdaMax = 10, prior = prior, priorType = priorType, interactionFlag = interactionFlag, precision=precision, precisionRatio = precisionRatio, priorInfo = priorInfo);
    BICResult[i] <- BIC(sample, result);
  } 
  return(BICResult);
}

BIC <- function(sample, result){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample); 
  gamma <- 0.5;
  likelihood <- n/2*(log(det(result$icov)) - sum(diag(sample$var %*% result$icov)) );
  edge <- (sum(result$adj) - p)/2*log(n);
  return(-2*likelihood + edge );
}