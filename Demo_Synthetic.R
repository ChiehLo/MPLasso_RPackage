run = TRUE;
eva = TRUE;



if (run ==TRUE){
  # Simulation Configuration
  OTUnum <- 50; # number of OTU
  simulations <- 2e2; # number of sample
  strength <- 0.2; # correlation strength
  numIt <- 2; # number of runs
  sampleLabel <- matrix(0, numIt, OTUnum*OTUnum);


  #MPLasso Configuration
  prior <- 0.5;
  precision <- 0.5;
  rmsePL = matrix(0, numIt);
  predPL = matrix(0, numIt, OTUnum*OTUnum);

  for (i in 1:numIt){

    #sample generation
    sample <- graph_select(OTUnum = OTUnum, Simulations = simulations, Strength = strength, type = "random");
    sampleLabel[i,] <- abs(as.vector(sample$adj));

    # model selection
    priorInfo <- generate_prior(sample, prior, precision)
    selectLambda <- selection(sample, type="MPLasso", prior = prior, precision = T, precisionRatio = precision, priorInfo = priorInfo);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;

    # optimal model
    resultPL <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", precision = T, precisionRatio = precision, priorInfo = priorInfo);
    predPL[i,] <- abs(as.vector(resultPL$icov));

    # calculate L1 distance
    evaluationPL <- evaluation(sample, resultPL);
    rmsePL[i] <- evaluationPL$RMSE_l0;
  }
}

if(eva==TRUE){
  tempLabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  tempPred <- split(t(predPL), rep(1:nrow(predPL), each = ncol(predPL)))
  rocPL <- ROC_new(tempLabel, tempPred)
  perfPL <- rocPL$perf;
  aucPL <- rocPL$auc;
  aucstdPL <- rocPL$aucsd
  acc <- acc_eval(tempLabel, tempPred);
  auprPL <- rocPL$aucpr
  auprstdPL <- rocPL$aucprsd
  print("MPLasso Evaluation Results: rmse (std), acc (std), aupr (std)")
  print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(rmsePL), sd(rmsePL), mean(acc), sd(acc), auprPL, auprstdPL))
}



