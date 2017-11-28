#####################################################################
# Real data preprocessing
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

read_OTU <- function(level = "genus", countPath, association = NULL, interaction = NULL){
  if(level == "genus"){
    rawData <- read.csv(countPath, sep = "\t", header=FALSE)
    rawData <- as.matrix(rawData)
    sample = count_data_process(rawData = t(rawData), association = association);
  }
  else if(level == "species"){
    rawData <- t(as.matrix(read.table(countPath)));
    sample = fraction_data_process(rawData, association, interaction);
  }
  return(sample);
}

count_data_process <- function(rawData, association){
  p <- ncol(rawData);
  data <- t(clr(rawData+0.1, 1));
  s <- var(data);
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(s[i,i]);
  }
  D <- diag(d);
  R <- D%*%s%*%D;
  A <- association;
  return(list("sample" = rawData, "fraction" = data , "var" = s, "cor" = R, "adj" = association));
}

fraction_data_process <- function(rawData, association, interaction){
  p <- ncol(rawData);
  data <- t(clr(rawData+0.001, 1));
  s <- var(data);
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(s[i,i]);
  }
  D <- diag(d);
  R <- D%*%s%*%D;
  A <- association;
  return(list("sample" = rawData, "fraction" = data , "var" = s, "cor" = R, "adj" = association, "inter" = interaction));
}
