#####################################################################
# Synthetic graph generation: random, AR4, hub, cluster, scale-free
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################



graph_select <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15, type = "random"){
  if(type=="random"){
    sample <- random_model(OTUnum, Simulations, Strength);
    #print("random");
  }
  else if(type == "AR4"){
    sample <- AR4_model(OTUnum, Simulations, Strength);
    #print("AR4");
  }
  else if(type == "hub"){
    sample <- hub_model(OTUnum, Simulations, Strength);
    #print("hub");
  }
  else if(type == "cluster"){
    sample <- cluster_model(OTUnum, Simulations, Strength);
    #print("cluster");
  }
  else if(type=="scale-free"){
    sample <- scale_free_model(OTUnum, Simulations, Strength);
    #print("scale-free");
  }
  return(sample);
}


AR4_model <- function(OTUnum = 50, Simulations = 1e2, MaxStrength = 0.15 ){
  sim <- huge::huge.generator(n = Simulations, d = OTUnum, graph = "band", v = 0.3, u = 0 , g= 4, verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

hub_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.2 ){
  sim <- huge::huge.generator(n = Simulations, d = OTUnum, graph = "hub", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

cluster_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  sim <- huge::huge.generator(n = Simulations, d = OTUnum, graph = "cluster", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

random_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  sim <- huge::huge.generator(n = Simulations, d = OTUnum, graph = "random", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

scale_free_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  sim <- huge::huge.generator(n = Simulations, d = OTUnum, graph = "scale-free", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}


