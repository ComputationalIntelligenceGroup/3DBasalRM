#' Learn BN
#'
#' Learn BN structure and parameters according to the training dataset. There are arcs from class to all the explanatory variables and cannot be in the other way direction.
#'
#' @param data_learn a dataset to learn the model
#' @param wl is TRUE if we force an arc from the class to all the features and FALSE otherwise
#' @param nump an integer. It is the maximum number of parents per node
#'
#' @return a list where its first element is the structure of the BN and the second is the parameters
#'
#' @export
learn_BN<-function(data_learn,wl = FALSE,nump=Inf)
{
  whitelist <- data.frame(class = "class",variables = colnames(data_learn))
  whitelist <- whitelist[-(which(whitelist[,2]=="class")),]

  blacklist <- data.frame(variables = colnames(data_learn), class = "class")
  blacklist <- blacklist[-(which(blacklist[,1]=="class")),]

  if (wl){
    bn <- bnlearn:::hc(data_learn, whitelist = whitelist, blacklist = blacklist, maxp = nump, restart = 1)
  } else{
    bn <- bnlearn:::hc(data_learn, blacklist = blacklist, maxp = nump, restart = 1)
  }
  parameters <- bnlearn:::bn.fit(bn, data_learn)

  return(list(structure = bn, params = parameters))
}


#' Predict with BN
#'
#' Predict the class of a dataset given a BN
#'
#' @param BN is an object of bnlearn containing the BN structure and parameters
#' @param pred_data a dataset to predict the classes of the instances
#'
#' @return a data.frame whose columns are the truth class, the probability of each class and the predicted class
#'
#' @export
pred_BN <- function(BN, pred_data)
{
  classes <- as.numeric(attributes(BN$params$class$prob)$dimnames[[1]])
  log_lik <- matrix(0, ncol = length(classes), nrow = nrow(pred_data))
  weights <- matrix(0, ncol = length(classes), nrow = nrow(pred_data))

  for(i in classes)
  {
    temp <- pred_data
    temp$class <- factor(i, levels = classes)
    log_lik[,which(i==classes)] <- bnlearn:::logLik.bn.fit(BN$params, temp,by.sample = T)
  }

  denominator <- logSumExp(log_lik)
  weights <- exp(sweep(log_lik, 1, denominator, "-"))
  pred_class<-classes[unlist(apply(log_lik,1,which.max))]

  #If the log-likelihood is under -50 the node is consider an outlier and it is denoted as a terminal node
  idx_outlier <- which(apply(log_lik, 1, max) < -50)

  #Save outliers as terminal
  if(length(idx_outlier) > 0)
  {
    for(idx in idx_outlier) {weights[idx,] <- c(1,0,0)}
  }

  #If there is a variable named class return the prediction
  if("class" %in% colnames(pred_data))
  {
    return(data.frame(truth = pred_data$class, prob = weights, resp = pred_class))
  }else{
    return(data.frame(prob = weights, resp = pred_class))
  }
}

#' Compute log sum exp function
#'
#' Compute log sum exp function to avoid underflow
#'
#' @param x is a matrix of log likelihood where each column is the log-likelihood of a class
#'
#' @return a vector with the probability of the observations
logSumExp<-function(x)
{
  y = apply(x, 1, max)
  x = sweep(x, 1, y, "-")
  s = y + log(rowSums(exp(x)))
  return(s)
}
