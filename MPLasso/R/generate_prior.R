generate_prior <- function(sample, prior, precision){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);
  adj <- sample$adj;
  rhoInter <- matrix(1, p, p);
  diag(rhoInter) <- 0
  rhoOccur <- matrix(0, p, p);
  # type 1: know the interaction, i.e., know adj[i,j] == 1 
  if (prior >= 0 && precision > 0){
    precisionRatio <- precision
    # 1 for 
    rhoInter <- matrix(1, p, p); #gLasso regularizer
    diag(rhoInter) <- 0
    # generate the constriant matrix 
    countPenalized <- 0
    totalCount <- 0
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 1 && runif(1, min=0, max=1) > prior && i != j){
          rhoInter[i,j] <- 0;
          if(runif(1, min=0, max=1) > (1- precisionRatio) ){
            countPenalized <- countPenalized + 1
            rhoInter[i,j] <- 1;
          }
          totalCount <- totalCount + 1
        }
      }
    }
    #print(countPenalized)
    #print(totalCount)
    while (countPenalized > 0){
      row = round(runif(1, min = 1, max = p))
      col = round(runif(1, min = 1, max = p))
      if(row!=col && rhoInter[row, col]==1){
        rhoInter[row,col] = 0;
        countPenalized <- countPenalized -1
      }
    }
    print(countPenalized)
  }
  
  # type 2: know non interacting edges, i.e., adj[i,j] == 0 (i.e., from co-occurrence method)
  if (prior >= 0 && precision > 0){
    precisionRatio <- precision
    rhoOccur <- matrix(0, p, p); #gLasso regularizer
    # generate the constriant matrix
    countPenalized <- 0
    totalCount <- 0
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior && i != j){
          rhoOccur[i,j] <- 1;
          if(runif(1, min=0, max=1) > (1- precisionRatio) ){
            countPenalized <- countPenalized + 1
            rhoOccur[i,j] <- 0;
          }
          totalCount <- totalCount + 1
        }
      }
    }
    print('type2:')
    print(totalCount)
    print(countPenalized)
    while(countPenalized > 0){
      row = round(runif(1, min = 1, max = p))
      col = round(runif(1, min = 1, max = p))
      if(row!=col && rhoOccur[row, col]==0){
        rhoOccur[row,col] = 1;
        countPenalized <- countPenalized -1
      }
    }
    print(countPenalized)
  }
  
  if(prior>=0 && precision == 0){
    rhoInter <- matrix(1, p, p); #gLasso regularizer
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 1 && runif(1, min=0, max=1) > prior && i!=j){
          rhoInter[i,j] <- 0;
        }
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior && i != j){
          rhoOccur[i,j] <- 1;
        }
      }
    }
  }
  return(list("rhoInter" = rhoInter, "rhoOccur" = rhoOccur));
}