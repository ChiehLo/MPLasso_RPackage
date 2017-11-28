# five body sites considered in this demo
bodySites <- c('Anterior_nares', 'Buccal_mucosa', 'Stool', 'Supragingival_plaque', 'Tongue_dorsum')
bodySites <- c('Anterior_nares')
# flag for evaluating different datasets
HMASM = F;
HMMCP = F;
HMQCP = T;

# reproducibility evaluation
repro = T;

if (HMASM == TRUE){
  for (i in 1:length(bodySites)){
    print("HMASM:")
    print(bodySites[i]);
    bodySite = bodySites[i];
    interactionFile = "/output/species_interaction.csv";
    interactionPath = paste("HMASM/", bodySite, interactionFile, sep= '');
    interactionPath <- system.file("extdata", interactionPath, package = "MPLasso")
    interactionMatrix <- interaction(interactionPath); # read interaction file

    occuPath = paste("HMASM/", bodySite, "/output/association.csv", sep = '');
    occuPath <- system.file("extdata", occuPath, package = "MPLasso")
    NoAssociation <- occurance(occuPath); # read occurance file

    countPath = paste("HMASM/", bodySite, "/output/SpeciesCount.txt", sep = '');
    countPath <- system.file("extdata", countPath, package = "MPLasso")
    sample <- read_OTU(level = "species", countPath = countPath, 1-NoAssociation, interactionMatrix); # read OTU

    prior = 0;
    interactionFlag = TRUE;
    selectLambda <- selection(sample, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
    resultPL <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);

    nInt <- 2;
    reproduceErrorArrayPL = matrix(0,nInt, 4);
    if (repro == TRUE){
      for (k in 1:nInt){
        numSample <- nrow(sample$sample)
        subSample <- sample;
        sampleIndex <- sort(sample.int(numSample, round(0.5*numSample)))
        subSample$sample <- subSample$sample[sampleIndex,]
        lassoSub <- fraction_data_process(subSample$sample, 1 - NoAssociation, interactionMatrix)
        selectLambda <- selection(lassoSub, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
        lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
        resultPLSample <- algorithm_select(lassoSub, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
        reproduceErrorPL <- repro_eval(resultPL$adj, resultPLSample$adj)
        reproduceErrorArrayPL[k,] <- reproduceErrorPL;

      }
      meanPL <- colMeans(1-reproduceErrorArrayPL)
      stdPL <- matrixStats::colSds(1-reproduceErrorArrayPL)
      print("MPLasso HMASM Evaluation Results: 25% (std), 50% (std), 75% (std), 100% (std)")
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", meanPL[1], stdPL[1], meanPL[2], stdPL[2], meanPL[3], stdPL[3], meanPL[4], stdPL[4]))
    }
  }
}

if (HMMCP == TRUE){
  for (i in 1:length(bodySites)){
    print("HMMCP:")
    print(bodySites[i]);
    bodySite = bodySites[i];
    occuPath = paste("HMMCP/v35/", bodySite, "/association.csv", sep = '');
    occuPath <- system.file("extdata", occuPath, package = "MPLasso")
    NoAssociation <- occurance(occuPath); # read occurance file


    countPath = paste("HMMCP/v35/", bodySite, "/genus_count_filter.csv", sep = '');
    countPath <- system.file("extdata", countPath, package = "MPLasso")
    sample <- read_OTU(level = "genus", countPath = countPath, 1 - NoAssociation);

    prior = 0;
    interactionFlag = FALSE;
    selectLambda <- selection(sample, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
    resultPLMM <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
    nInt <- 2;
    reproduceErrorArrayPL = matrix(0,nInt, 4);
    if (repro == TRUE){
      for (k in 1:nInt){
        numSample <- nrow(sample$sample)
        subSample <- sample;
        sampleIndex <- sort(sample.int(numSample, round(0.5*numSample)))
        subSample$sample <- subSample$sample[sampleIndex,]
        lassoSub <- count_data_process(subSample$sample, 1 - NoAssociation)
        selectLambda <- selection(lassoSub, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
        lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
        resultPLSample <- algorithm_select(lassoSub, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
        reproduceErrorPL <- repro_eval(resultPLMM$adj, resultPLSample$adj)
        reproduceErrorArrayPL[k,] <- reproduceErrorPL;
      }
      meanPL <- colMeans(1-reproduceErrorArrayPL)
      stdPL <- matrixStats::colSds(1-reproduceErrorArrayPL)
      print("MPLasso HMMCP Evaluation Results: 25% (std), 50% (std), 75% (std), 100% (std)")
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", meanPL[1], stdPL[1], meanPL[2], stdPL[2], meanPL[3], stdPL[3], meanPL[4], stdPL[4]))
    }

  }
}

if (HMQCP == TRUE){
  for (i in 1:length(bodySites)){
    print("HMQCP:")
    print(bodySites[i]);
    bodySite = bodySites[i];
    occuPath = paste("HMQCP/v35/", bodySite, "/association.csv", sep = '');
    occuPath <- system.file("extdata", occuPath, package = "MPLasso")
    noAssociation <- occurance(occuPath);

    countPath = paste("HMMCP/v35/", bodySite, "/genus_count_filter.csv", sep = '');
    countPath <- system.file("extdata", countPath, package = "MPLasso")
    sample <- read_OTU(level = "genus", countPath = countPath, association = 1 - noAssociation);

    prior = 0;
    interactionFlag = FALSE;
    selectLambda <- selection(sample, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
    resultPLMQ <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
    nInt <- 2;
    reproduceErrorArrayPL = matrix(0,nInt, 4);
    if (repro == TRUE){
      for (k in 1:nInt){
        numSample <- nrow(sample$sample)
        subSample <- sample;
        sampleIndex <- sort(sample.int(numSample, round(0.5*numSample)))
        subSample$sample <- subSample$sample[sampleIndex,]
        lassoSub <- count_data_process(subSample$sample, 1 - noAssociation)
        selectLambda <- selection(lassoSub, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
        lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
        resultPLSample <- algorithm_select(lassoSub, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
        reproduceErrorPL <- repro_eval(resultPLMQ$adj, resultPLSample$adj)
        reproduceErrorArrayPL[k,] <- reproduceErrorPL;

      }
      meanPL <- colMeans(1-reproduceErrorArrayPL)
      stdPL <- matrixStats::colSds(1-reproduceErrorArrayPL)
      print("MPLasso HMQCP Evaluation Results: 25% (std), 50% (std), 75% (std), 100% (std)")
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", meanPL[1], stdPL[1], meanPL[2], stdPL[2], meanPL[3], stdPL[3], meanPL[4], stdPL[4]))
    }

  }
}
