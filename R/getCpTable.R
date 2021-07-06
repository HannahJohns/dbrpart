#' Get Complexity Parameter Table for a dbrpart fit
#'
#' @param x an object of class \code{"Dbrpart"}
#'
#' @returns a list containing
#' \itemize{
#'     \item \code{cpTable}
#'     \item \code{index} - the row of cpTable corresponding to the best fit using the 1-SE rule 
#' }
#'
#' @export
getCpTable <- function(x)
{
  if(x$controls$stopMethod!="cv") stop("Must use crossvalidation for CpTable to be meaningful!")
  
  loss <- x$loss

  stopifnot(length(x$b)==ncol(x$cv_prediction))
  
  as.data.frame(t(sapply(1:length(x$b),function(i){
    
    beta <- x$b[i]
    nTerminal <- length(unique(x$classTree[,i]))
    
    sapply(1:length(x$y),function(j){
      loss[
        as.numeric(x$y)[j],
        x$prediction[j,i]+1
      ]  
    }) -> IndividualLoss
    
    MeanLoss <- sum(IndividualLoss)/length(x$y)
    
    sapply(1:length(x$y),function(j){
      loss[
        as.numeric(x$y)[j],
        x$cv_prediction[j,i]+1
      ]  
    }) -> xIndividualLoss
    
    xMeanLoss <- sum(xIndividualLoss)/length(x$y)
    
    varLoss <- sum((xIndividualLoss - xMeanLoss)^2)/length(x$y)
    xSE <- sqrt(varLoss/length(x$y))
    
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
