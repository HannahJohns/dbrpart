#' Control for dbrpart fits
#'
#'
#' @param minSplit the minimum number of observations that must exist in a node in order for a split to be attempted.
#' 
#' @param minBucket the minimum number of observations in any terminal node
#' 
#' @param classMethod the method for classifying members in a node
#' 
#' @param maxDepth the maximum depth of the tree
#' 
#' @param distanceAdjust a penalty when deciding bubble center based on total distance to other observations.
#' 
#' @param trace the level of information to print while dbrpart runs
#' 
#' @param stopMethod the method for determining when stopping occurs
#' 
#' @param nSims Number of permutations to perform in permutation stopping criterion
#' 
#' @param signifLevel Significance level for permutation stopping criterion
#' 
#' @param nCenters Number of bubble centers to test
#'
#' @param ... 
#'
#' @export
dbrpart.control <- function(minSplit = 20,
                            minBucket = round(minSplit/3),
                            classMethod="default", # c("mean","median","mode"),
                            maxDepth = 20,
                            distanceAdjust=0,
                            trace = 3,
                            stopMethod = 0, #c("cv","permutation","exhaust")
                            nSims=5000,
                            signifLevel=0.05,
                            pCenters=1,
                            combinationPower = 2
                            )
{
  
  if(2*minBucket > minSplit)
  {
    warning("minBucket should be less than minSplit/2")
  }
  
  
  
  out <- list(
    minSplit=minSplit,
    minBucket=minBucket,
    classMethod=classMethod,
    maxDepth=maxDepth,
    distanceAdjust=distanceAdjust,
    trace=trace,
    stopMethod=stopMethod,
    nSims=nSims,
    signifLevel=signifLevel,
    pCenters=pCenters,
    combinationPower = combinationPower
  )
  
  class(out) <- c(class(out),"dbrpart.control")
  return(out)
  
}