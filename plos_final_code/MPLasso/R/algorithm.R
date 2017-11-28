#####################################################################
# MPLasso implementation
# You can implement your algorithm here
# @author Chieh Lo
# @date 05/13/2017
#####################################################################



algorithm_select <- function(sample, lambdaMin = 0.01, lambdaMax = 10, prior = 1, type = "MPLasso", priorType = TRUE, interactionFlag = FALSE, precision = FALSE, precisionRatio, priorInfo = priorInfo){
  if(type=="MPLasso"){
    result <- MPLasso(sample, lambdaMin, lambdaMax, prior = prior, priorType = priorType, interactionFlag = interactionFlag, precision = precision, precisionRatio = precisionRatio, priorInfo = priorInfo)
  }
  return(result);
}


MPLasso <- function(sample, lambdaMin = 0.01, lambdaMax = 10, prior = 1, priorType, interactionFlag, precision = FALSE, precisionRatio = precisionRatio, priorInfo = priorInfo){
  s <- sample$var;
  adj <- sample$adj;
  if(interactionFlag==TRUE){
    interaction <- sample$inter;
  }
  p <- nrow(s);
  rho <- lambdaMin*matrix(1, p, p); #gLasso regularizer

  ## no pricision ratio (i.e., all priors are correct)
  ## adj encodes true graph structure. For synthetic data, we use the true graph structure as our prior info.
  ## For real data, adj encodes the co-occurrence information (i.e., taxa is not associated with each other)
  if (priorType == TRUE && precision == FALSE){
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior && i>j){
          rho[i,j] <- lambdaMax;
          rho[j,i] <- lambdaMax;
        }
      }
    }
  }

  ## consider pricision ratio (given prior information)
  if (priorType == TRUE && precision == TRUE && prior < 1){
    rhoInter <- priorInfo$rhoInter
    rhoOccur <- priorInfo$rhoOccur
    rho <- lambdaMin*matrix(1, p, p); #gLasso regularizer
    if (precisionRatio==0){
      lambdaMax = 1
    }
    else{
      lambdaMax = lambdaMin/precisionRatio
    }

    for (i in 1:p){
      for (j in 1:p){
        if(rhoOccur[i,j]==1 && rhoInter[i,j]==1){
          rho[i,j] <- lambdaMax
          rho[j,i] <- lambdaMax
        }
      }
    }
  }
  if(prior == 1){
    rho <- lambda_min*matrix(1, p, p); #gLasso regularizer
  }




  if(interactionFlag == TRUE){
    for (i in 1:p){
      for (j in 1:p){
        if (abs(interaction[i,j]) == 1 && runif(1, min=0, max=1) > prior && i!=j){
          rho[i,j] <- 0.01;
          rho[j,i] <- 0.01;
        }
      }
    }
  }
  a<-glasso::glasso(s, rho=rho)
  covOpt <- a$w;
  icovOpt <- a$wi;
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(covOpt[i,i]);
  }
  D <- diag(d);
  corrOpt <- D%*%covOpt%*%D;

  adjOpt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(icovOpt[i,j]) > 0){
        adjOpt[i,j] <- 1;
        adjOpt[j,i] <- 1;
      }
    }
  }


  return(list("cor" = corrOpt, "cov" = covOpt, "icov" = icovOpt, "adj" = adjOpt ));
}



