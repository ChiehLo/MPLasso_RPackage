#####################################################################
# Performance evaluation
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

evaluation <- function(sample, result){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);

  RMSE <- RMSE_l0(sample$cor, result$cor);
  ICOV <- ICOV_l0(MASS::ginv(sample$cor), result$icov);
  RMSE_f <- RMSE_F(sample$cor, result$cor);
  ICOV_f <- ICOV_F(MASS::ginv(sample$cor), result$icov);

  return(list("RMSE_l0" = RMSE, "ICOV_l0" = ICOV, "RMSE_F" = RMSE_f, "ICOV_F" = ICOV_f));
}

RMSE_l0 <- function(true_cor, estimate_cor){
  p <- nrow(true_cor);
  Error <- abs(true_cor - estimate_cor);
  RMSE <- sum(sum(Error))/(p*(p-1));
  return(RMSE);
}

ICOV_l0 <- function(true_icov, estimate_icov){
  p <- nrow(true_icov);
  Error <- abs(true_icov - estimate_icov);
  RMSE <- sum(sum(Error))/(p*(p-1));
  return(RMSE);
}


RMSE_F <- function(true_cor, estimate_cor){
  p <- nrow(true_cor);
  Error <- abs(true_cor - estimate_cor)*abs(true_cor - estimate_cor);
  RMSE <- sqrt(sum(sum(Error)));
  return(RMSE);
}

ICOV_F <- function(true_icov, estimate_icov){
  p <- nrow(true_icov);
  Error <- abs(true_icov - estimate_icov)*abs(true_icov - estimate_icov);
  RMSE <- sqrt(sum(sum(Error)));
  return(RMSE);
}


ROC_eval <- function(label, pred){

  pred <- prediction(pred, label)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred, 'auc')
  aucstd <- sd(simplify2array(auc@y.values))
  auc <- mean(simplify2array(auc@y.values))
  return(list("perf" = perf, "auc" = auc, "aucsd" = aucstd));
}

ROC_new <- function(label, pred){

  pred <- ROCR::prediction(pred, label)
  perf <- ROCR::performance(pred,"tpr","fpr")
  auc <- ROCR::performance(pred, 'auc')
  aucstd <- sd(simplify2array(auc@y.values))
  auc <- mean(simplify2array(auc@y.values))
  f <- ROCR::performance(pred,"prec","rec")
  x <- f@x.values
  y <- f@y.values

  AUCPr = matrix(0,length(x), 1);
  for (i in 1:length(x)){
    xTemp <- as.matrix(x[[i]])
    yTemp <- as.matrix(y[[i]])
    yTemp[1] <- yTemp[2];
    AUCPr[i] = pracma::trapz(xTemp, yTemp)
    f@y.values[[i]][1] <- f@y.values[[i]][2];
  }
  aucpr <- mean(AUCPr)
  aucprsd <- sd(AUCPr)

  return(list("perf" = perf, "auc" = auc, "aucsd" = aucstd, "f" = f, "aucpr" = aucpr, "aucprsd" = aucprsd));
}


acc_eval <- function(label, pred){
  pred <- ROCR::prediction(pred, label)
  perf <- ROCR::performance(pred,"acc")
  acc <- matrix(0, length(label))
  for (i in 1:length(label)){
    index <- which.max(perf@y.values[[i]])
    acc[i] <- perf@y.values[[i]][index];
  }
  return(acc);
}

repro_eval <- function(true_adj, test_adj){
  row_sum <- rowSums(true_adj);
  degree <- sort(row_sum, index.return = TRUE, decreasing=TRUE);
  degree_index <- degree$ix;
  num_taxa <- nrow(true_adj);
  percent = c(0.25, 0.5, 0.75, 1);
  reproError = matrix(0, 1, length(percent))
  for(i in 1:length(percent)){
    temp_index <- degree_index[1:round(num_taxa*percent[i])]
    #print(temp_index)
    num_edge_true <- sum(true_adj[temp_index,]) - round(num_taxa*percent[i])
    num_edge_test <- sum(test_adj[temp_index,]) - round(num_taxa*percent[i])
    error <- true_adj[temp_index,] -test_adj[temp_index,];
    #print(nrow(error))
    reproError[i] <- sum(abs(error))/round((num_taxa*num_taxa*percent[i]))
    #print(sum(abs(error)))
    #print(round(num_taxa*percent[i]))
  }

  return(reproError);

}
