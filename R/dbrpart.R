#' Recursive composite k-nearest neighbour partitioning
#'
#' Variation on CART using distance-based paritioning
#'
#' @param target Response Variable. Must be either a numeric, logical or factor vector, or a distance matrix.
#'
#' @param distanceList A list of distance matrices. rcknnpart will partition using these features.
#'
#' @param priors Prior distribution for categorical outcome frequency. Ignored if target is not a factor.
#' 
#' @param loss Loss matrix for categorical outcome frequency. Ignored if target is not a factor.
#' @param factorList 
#' @param numericList 
#' @param distanceCombinations 
#' @param alteredPriors 
#' @param controls 
#' @param soma.controls 
#' 
#' @returns A list containing the following
#'
#' \itemize{
#' 
#' If cross-validation is being used:
#' 
#' \item `a` - the alpha complexity thresholds for each prune
#' \item `b` - beta, the geometric mean between successive pruning values
#' \item `nTerm` - the number of terminal nodes at pruning step
#' \item `model_risk` - the average error at each pruning step
#' \item `risk` - the sum of errors at each pruning step
#' \item `cv_risk` - the sum of errors at each pruning step under cross-validation
#' \item `classTree` - the terminal node that each observation belongs to at each pruning step
#' \item `prediction` - the predicted value of each observation at each pruning step (under the model)
#' \item `details`
#' \item `cv_prediction` - the predicted value of each observation at each pruning step (under cross-validation)
#' \item `cv_details`
#' 
#' For all methods, the following:
#' 
#' \item `tree` -  a list containing the details of each node in the tree:
#' \itemize{
#'     \item `label` - the label for the node. This is currently zero indexed and should be changed.
#'     \item `parent` - the label for the node's parent. Currently, this doesn't contain anything and should probably be removed.
#'     \item `probability` - the probability of a random patient being in this node
#'     \item `groupMembers` - a list of observations in this node. This is currently zero indexed and should be changed.
#'     \item `groupSize` - the count of observations in each node.
#'     \item `risk` - the ``risk'' of the node. This changes depending on the target:
#'     \itemize{
#'         \item If the target is nominal or ordinal, then it contains the 
#'         \item If the target is scalar, then it should contain 
#'         \item If the target is distance-based,, then
#'     }
#'     \item `pruneAtStep`
#'     \item `impurity`
#'     \item `classif`
#'     \item `mean`
#'     \item `terminal`
#'     \item `alpha`
#'     \item `splitMethod`
#'     \item `bestSplitIndex`
#'     \item `splitName`
#'     \item `splitUtility`
#'     \item `goodnessOfFit`
#'     \item `center`
#'     \item `radius`
#'     \item `pVal`
#'     \item `inChild`
#'     \item `outChild`
#' }
#' 
#' \item `y` - the target variable.
#' \item `featureNames`
#' \item `classNames`
#' \item `mode`
#' \item `controls`
#' \item `foldList`
#' \item `bestPruneStep`
#' } 
#'
#' 
#' 
#' 
#' @export
dbrpart <- function(target,
                    factorList = NULL,
                    numericList=NULL,
                    distanceList=NULL,
                    distanceCombinations=NULL,
                    priors = NULL,
                    alteredPriors = FALSE,
                    loss = NULL,
                    controls=NULL,
                    soma.controls=NULL
                    )
{
  
  
  
  
  # Before anything else, determine if we're trying to build:
  if("factor" %in% class(target) | "logical" %in% class(target))
  { # A Classification Tree

    # Code response for Rcpp
    target_discrete <- factor(target)
    target_continuous <- 0
    target_distance <- matrix(0)
    mode <- 0

    # Swap to zero indexing
    target_discrete <- as.numeric(target_discrete)-1
    nObs <- length(target_discrete)

    # Code metadata about response for Rcpp
    classNames <- levels(target)
    nCategories <- length(classNames)

    if(is.null(priors))
    {
      classFreq <- as.vector(table(target))
      priors<- classFreq/sum(classFreq)
    }

    if(is.null(loss))
    {
      if("ordered" %in% class(target))
      {
        # If factor is ordered, use linear weighting
        loss <- diag(nCategories)
        loss <- abs(row(loss)-col(loss))
      }
      else
      {
        # Otherwise, use constant misclassification cost
        loss <- 1-diag(nCategories)
      }
    }

  }
  else if (class(target) == "numeric")
  { # A Regression Tree
    target_discrete <- 0
    target_continuous <- as.numeric(target)
    nObs <- length(target_continuous)
    
    target_distance <- matrix(0)
    mode <- 1
  }
  else if (class(target) == "dist")
  { # A Distance-Based Regression Tree
    target_discrete <- 0
    target_continuous <- 0
    target_distance <- as.matrix(target)
    nObs <- length(target_distance)
    
    
    mode <- 2
  }
  else
  { # Some sort of unsupported nonsense
    stop("target is not a recognised class! Must be factor, logical, numeric or dist")
  }

  
  if(!is.null(factorList))
  {
    if(prod(sapply(factorList, function(x){"factor" %in% class(x)})) != 1)
    {
      stop("factorList is not a list of class factor")
    }
  }
  
  
  if(!is.null(numericList))
  {
    if(prod(sapply(numericList, function(x){"numeric" %in% class(x)})) != 1)
    {
      stop("numericList is not a list of class numeric")
    }
  }
  
  
  if(!is.null(distanceList))
  {
    if(prod(sapply(distanceList, function(x){"dist" %in% class(x)})) != 1)
    {
      stop("distanceList is not a list of class dist")
    }
    
    if(is.null(distanceCombinations))
    {
      distanceCombinations <- as.list( (1:length(distanceList))-1 )
    }
  }
  
  
  if(!is.null(priors))
  {
    if(mode!=0)
    {
      warning("priors option ignored where target is not a factor")
    }
  }
  
  
  if(alteredPriors)
  {
    if(mode!=0)
    {
      warning("alteredpriors ignored where target is not a factor")
    }
  }
  
  
  if(!is.null(loss))
  {
    if(mode!=0)
    {
      warning("loss ignored where target is not a factor")
    }
  }
  
  
  if(!is.null(controls))
  {
    if(!("dbrpart.control" %in% class(controls)))
    {
      stop("controls should be specified by function dbrpart.control()")
    }
  }
  
  
  # Check outputs to see what we're using as a basis for tree selection
  if(is.null(distanceList) & is.null(numericList) & is.null(factorList))
  {
    stop("Must provide at least one of distanceList, numericList or factorList")
  }
  
  
  if(is.null(distanceList))
  {
    distanceList <- list()
  }
  if(is.null(numericList))
  {
    numericList <- list()
  }
  if(is.null(factorList))
  {
    factorList <- list()
  }
  

  
  # Now that we know what mode we're running in, we'll grab parameters from dbrpart.control
  if(is.null(controls))
  {
    controls <- dbrpart.control()
  }
  
  # Take stopMethod and give it a magic number to pass to c++
  if(controls$stopMethod=="default")
  {
    # The default method will depend on what we're trying to fit the tree on
    controls$stopMethod <- c("cv","cv","permutation")[mode+1]
  }
  
  
  if(controls$classMethod=="mean")
  {
    if(mode == 0)
    {
      stop("classMethod of mean doesn't make sense for factor target")
    }
    else if (mode==2)
    {
      warning("classMethod ignored for distance target")
    }
    
    classMethod <- 0
  }
  else if(controls$classMethod=="median")
  {
    if(mode==0 & !("ordered" %in% class(target)))
    {
      stop("classMethod of median does not make sense if target is not an ordered factor")
    }
    else if(mode == 1)
    {
      stop("classMethod of median not implemented for numeric targets")
    }
    else if (mode==2)
    {
      warning("classMethod ignored for distance target")
    }
    classMethod <- 1
  }
  else if(controls$classMethod=="mode")
  {
    if(mode == 1)
    {
      stop("classMethod of median not implemented for numeric targets")
    }
    else if (mode==2)
    {
      warning("classMethod ignored for distance target")
    }
    
    classMethod <- 2
  }
  else if(controls$classMethod=="default")
  {
    if(mode==0)
    {
      if("ordered" %in% class(target))
      {
        classMethod <- 1 # Median
      }
      else
      {
        classMethod <- 2 # Mode
      }
    }
    if(mode == 1)
    {
      classMethod <- 0 # Mean
    }
    else if (mode==2)
    {
      classMethod <- 0 # Mean (dummy)
    }
  }
  
  if(is.null(soma.controls))
  {
    soma.controls <- soma.control()
  }
  
  
  minSplit <- controls$minSplit
  minBucket <- controls$minBucket
  maxDepth <- controls$maxDepth
  trace <- controls$trace
  stopMethod <- stopMethod <- c(cv=0,permutation=1,none=2)[controls$stopMethod]
  nSims <- controls$nSims
  signifLevel <- controls$signifLevel
  pCenters <- controls$pCenters
  combinationPower <- controls$combinationPower
  
  
  soma_pathLength <- soma.controls$pathLength
  soma_step <- soma.controls$step
  soma_popSize <- soma.controls$popSize
  soma_nDirections <- soma.controls$nDirections
  soma_migrations <- soma.controls$migrations
  soma_minDiv <- soma.controls$minDiv
  soma_PRT <- soma.controls$PRT
  soma_seed <- soma.controls$seed
  
  
   
  if(mode!=0)
  {
    # We don't need classification tree features,
    # but need to supply something.
    classNames <- 0
    nCategories <- 0
    priors<- 0
    loss <- matrix(0)
  }
  
  # With response variable sorted, check and clean independent variables
  featureNames <- names(distanceList)
  featureIDs <- rownames(distanceList[[1]])

  distanceList <- lapply(distanceList,function(x){as.matrix(x)})
  testNA <- sapply(distanceList,function(x){sum(is.na(x))>0})
  if(sum(testNA)>0)
  {
    stop(sprintf("NA in %s",paste(names(testNA)[testNA],collapse=",")))
  }
  
  # Get half of minimum nonzero distance for each feature
  # Used to prevent shennanigans with floating point comparisons
  acceptableError <- sapply(distanceList,function(x){min(x[x>0],na.rm = TRUE)})/2

  # Generate fold list for cross-validation.
  # Generated at the R level so that set.seed() impacts it.
  foldList <- rep((1:10)-1, length.out=length(target))
  foldList <- foldList[order(runif(length(foldList)))]

  fit_out <- fit_dbrpart(in_target_discrete=target_discrete,
                     in_target_continuous=target_continuous,
                     in_target_distance=target_distance,
                     in_factorList=factorList,
                     in_numericList=numericList,
                     in_distanceList=distanceList,
                     in_distanceCombinations=distanceCombinations,
                     in_combinationPower = combinationPower,
                     in_acceptable_error=acceptableError,
                     in_mode=mode,
                     in_stopMethod=stopMethod,
                     in_nSims=nSims,
                     in_signifLevel=signifLevel,
                     in_fold=foldList,
                     in_nCategories = nCategories,
                     in_priors = priors,
                     in_loss = loss,
                     in_minSplit = minSplit,
                     in_minBucket = minBucket,
                     in_maxDepth = maxDepth,
                     in_alteredPriors=alteredPriors,
                     in_trace = trace,
                     in_classMethod = classMethod,
                     in_pCenters = pCenters,
                     soma_pathLength=soma_pathLength,
                     soma_step=soma_step,
                     soma_PRT=soma_PRT,
                     soma_popSize=soma_popSize,
                     soma_nDirections=soma_nDirections,
                     soma_migrations=soma_migrations,
                     soma_minDiv=soma_minDiv,
                     soma_seed = soma_seed)

  
  # We don't need to return everything here
  
  out <- fit_out
  
  
  out$y <- target
  out$featureNames <- featureNames
  out$featureIDs <-featureIDs
  out$classNames <- classNames
  
  # Magic number codes aren't particularly useful as outputs, interpret them
  out$mode <- c("classification","regression","db-regression")[mode+1]
  out$controls <- controls
  
  if(out$mode=="classification")
  {
    out$loss=loss
  }
  
  
  
  out$foldList <- foldList
  
  out$bestPruneStep <- length(out$b)-min(which(out$cv_risk==min(out$cv_risk)))
  
  class(out) <- c("Dbrpart",class(out))

  return(out)
}


#' Print a Dbrpart tree
#'
#' @param ... 
#' @param x
#' 
#'
#'@export
print.Dbrpart <- function(x,...)
{
  parms <- list(...)
  
  if(is.null(parms$pruneStep))
  {
    if(x$controls$stopMethod=="cv")
    {
      cpTable <- getCpTable(x)
      parms$pruneStep <- cpTable$cpTable[cpTable$index,"pruneStep"]
    }
    else
    {
      parms$pruneStep <- -Inf 
    }
  }
  
  if(is.null(parms$label))
  {
    parms$label <- 1
  }
  
  if(x$mode=="classification")
  {
    cat(sprintf("%s%d) n=%d, r=%1.2f, class=%d",
                paste(rep(" ",x$tree$depths[parms$label]),collapse=""),
                x$tree$label[parms$label],
                x$tree$groupSize[parms$label],
                x$tree$risk[parms$label],
                x$tree$classif[parms$label]
    ))
  }
  else if(x$mode=="db-regression")
  {
    cat(sprintf("%s%d) n=%d, r=%1.2f",
                paste(rep(" ",x$tree$depths[parms$label]),collapse=""),
                x$tree$label[parms$label],
                x$tree$groupSize[parms$label],
                x$tree$risk[parms$label]
    ))
  }
  
  if(x$tree$terminal[parms$label] | x$tree$pruneAtStep[parms$label] <= (parms$pruneStep))
  {
    cat(" *\n")
  }
  else
  {
    cat(sprintf("| s:%s c:%d r:%f",
                x$tree$splitName[parms$label],
                x$tree$center[parms$label],
                x$tree$radius[parms$label]
    ))
    
    if(x$controls$stopMethod=="permutation")
    {
      cat(sprintf(" (p=%0.4f)", x$tree$pVal[parms$label]))
    }
    
    cat("\n")
    print.Dbrpart(x,label=x$tree$inChild[parms$label]+1,pruneShift = parms$pruneShift )
    print.Dbrpart(x,label=x$tree$outChild[parms$label]+1, pruneShift = parms$pruneShift )
  }
}

