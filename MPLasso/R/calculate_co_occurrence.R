#####################################################################
# Calculate the prior information
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

occurance <- function(occuPath, totalPaper = 615268){
  a <- read.csv(occuPath, sep = "\t", header=FALSE)
  countTable <- as.matrix(a)
  n = nrow(countTable);
  pValue = matrix(0, n, n);
  CILow = matrix(0, n, n);
  CIHigh = matrix(0, n, n);
  totalPaper = sum(diag(countTable))
  for (i in 1:n){
    for (j in 1:n){
      if (j>i){
        contingency <- matrix(c(countTable[i,j], countTable[j,j] , countTable[i, i], totalPaper), nrow = 2)
        temp <- fisher.test(contingency, alternative = "two.sided");
        pValue[i,j] <- temp$p.value;
        CILow[i,j] <- temp$conf.int[[1]];
        CIHigh[i,j] <- temp$conf.int[[2]];
      }
    }
  }
  
  noAssociation <- matrix(0,n,n);
  adjust <- round(p.adjust(pValue, method = "bonferroni"), 5)
  corrected <- matrix(adjust, ncol=n, nrow = n);
  for (i in 1:n){
    for (j in 1:n){
      if (j>i && (CIHigh[i,j] > 1 && CILow[i,j] < 1) && corrected[i,j] == 1){
        noAssociation[i,j] <- 1;
        noAssociation[j,i] <- 1;
      }
    }
  }
  return(noAssociation)
}

interaction <- function(interactionPath){
  a <- read.csv(interactionPath, sep = "\t", header=FALSE)
  interaction <- as.matrix(a)
  return(interaction)
}
