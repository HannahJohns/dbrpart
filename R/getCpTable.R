
getCpTable <- function(fit)
{
  
  
  # Loss matrix wasn't stored oops lol
  loss <- diag(length(unique(fit$y)))
  loss <- abs(row(loss)-col(loss))
  
  stopifnot(length(fit$b)==ncol(fit$cv_prediction))
  
  as.data.frame(t(sapply(1:length(fit$b),function(i){
    
    beta <- fit$b[i]
    nTerminal <- length(unique(fit$classTree[,i]))
    
    sapply(1:length(fit$y),function(j){
      loss[
        as.numeric(fit$y)[j],
        fit$prediction[j,i]+1
      ]  
    }) -> IndividualLoss
    
    MeanLoss <- sum(IndividualLoss)/length(fit$y)
    
    sapply(1:length(fit$y),function(j){
      loss[
        as.numeric(fit$y)[j],
        fit$cv_prediction[j,i]+1
      ]  
    }) -> xIndividualLoss
    
    xMeanLoss <- sum(xIndividualLoss)/length(fit$y)
    
    varLoss <- sum((xIndividualLoss - xMeanLoss)^2)/length(fit$y)
    xSE <- sqrt(varLoss/length(fit$y))
    
    return(c(pruneStep = i-1,beta=beta, nTerm = nTerminal, err=MeanLoss,xerr=xMeanLoss,xse=xSE))
    
  }))) -> cpTable
  
  cpTable$xse <- cpTable$xse/max(cpTable$xerr)
  cpTable$xerr <- cpTable$xerr/max(cpTable$xerr)
  cpTable$err <- cpTable$err/max(cpTable$err)
  
  threshold <- cpTable[which(cpTable$xerr==min(cpTable$xerr)),]
  threshold <- threshold[which(threshold$nTerm==min(threshold$nTerm)),]
  threshold <- threshold$xerr+threshold$xse
  
  chosenTree <- cpTable[which(cpTable$xerr <= threshold),]
  chosenTree <- chosenTree[which(chosenTree$nTerm == min(chosenTree$nTerm)),]
  
  return(list(cpTable = cpTable, index=chosenTree$pruneStep+1))
  
}
