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
  stopMethod <- controls$stopMethod
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

  out <- fit_dbrpart(in_target_discrete=target_discrete,
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

  out$y <- target
  out$featureNames <- featureNames
  out$featureIDs <-featureIDs
  out$classNames <- classNames
  out$mode <- c("classification","regression","db-regression")[mode+1]
  
  out$foldList <- foldList
  
  out$bestPruneStep <- length(out$b)-min(which(out$cv_risk==min(out$cv_risk)))
  
  class(out) <- c("Dbrpart",class(out))

  return(out)
}



print.Dbrpart <- function(x,...)
{
  parms <- list(...)
  
  if(is.null(parms$pruneShift))
  {
    parms$pruneShift <- -Inf # This is currently broken lol
  }
  
  if(is.null(parms$label))
  {
    parms$label <- 1
  }
  
  cat(sprintf("%s%d) n=%d, r=%1.2f, class=%d",
              paste(rep(" ",x$tree$depths[parms$label]),collapse=""),
              x$tree$label[parms$label],
              x$tree$groupSize[parms$label],
              x$tree$risk[parms$label],
              x$tree$classif[parms$label]
  ))
  
  if(x$tree$terminal[parms$label] | x$tree$pruneAtStep[parms$label] <= (x$bestPruneStep+parms$pruneShift))
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
    
    cat("\n")
    print.Dbrpart(x,label=x$tree$inChild[parms$label]+1,pruneShift = parms$pruneShift )
    print.Dbrpart(x,label=x$tree$outChild[parms$label]+1, pruneShift = parms$pruneShift )
  }
}

# 
# print.Dbrpart <- function(x,...){
# 
#   parms <- list(...)
#   
#   
#   if(is.null(parms$label)) # Label of null indicates start of tree
#   {
#     parms$label <- 0
#   }
#   
#   if(is.null(parms$beta)) # If beta is not given, use best under CV
#   {
#     parms$beta <- min(x$b[x$cv_risk == min(x$cv_risk)])
#   }  
# 
#   
#   if(x$mode=="classification")
#   {
# 
#     cat(sprintf("   %s\n",
#                 paste(
#                       substring(
#                         paste(parms$classNames,"                 ",sep=""),
#                         0,10
#                       ),
#                       collapse=" ")
#                 )
#         )
# 
#     cat(sprintf("%s: %s | %s\n",
#                 x$tree$label,
#                 paste(sprintf("%0.8f",x$tree$classProb),collapse=" "),
#                 parms$featureNames[x$tree$bestD+1]
#                 )
#     )
#     
#   }
#   else if(parms$mode=="regression")
#   {
#     cat(sprintf("%s: %s | %s\n",
#                 x$tree$label,
#                 sprintf("%0.2f (%0.2f)",
#                           mean(parms$y[x$tree$groupMembers+1]),
#                           sd(parms$y[x$tree$groupMembers+1])
#                           ),
#                   parms$featureNames[x$tree$bestD+1]
#       ))
#   }
#   else
#   {
#     cat("uhhhh")
#   }
# 
# 
#     for(i in 1:length(x$tree$subtree))
#     {
#       cat(rep("  ",parms$depth))
#       print.Rcknnpart(x$tree$subtree[[i]],depth=parms$depth+1,
#                       beta=parms$beta,
#                       featureNames=parms$featureNames,
#                       mode=parms$mode,
#                       y=parms$y)
#     }
#   }
#   else # This function was called recursively
#   {
#     splitDecision <- "TERMINATE"
#     if(!is.null(x$bestD))
#     {
#       splitDecision <- parms$featureNames[x$bestD+1]
#     }
# 
#     cat(rep("  ",parms$depth))
# 
# 
#     if(parms$mode=="classification")
#     {
#       cat(sprintf("%s: %s | %s\n",
#                   x$label,
#                   paste(sprintf("%0.3f",x$classProb),collapse=" "),
#                   splitDecision
#       ))
#     }
#     else if(parms$mode=="regression")
#     {
#       cat(sprintf("%s: %s | %s\n",
#                   x$label,
#                   sprintf("%0.2f (%0.2f)",
#                           mean(parms$y[x$groupMembers+1]),
#                           sd(parms$y[x$groupMembers+1])
#                           ),
#                   splitDecision
#       ))
#     }
#     else
#     {
#       cat(sprintf("%s | %s\n",
#                   x$label,
#                   splitDecision
#       ))
#     }
# 
#     pruned <- ifelse(is.null(x$criticalAlpha),
#                      x$terminal,
#                      ((x$terminal) | (parms$beta>=x$criticalAlpha))
#     )
# 
#     if(!pruned)
#     {
#       cat(rep("  ",parms$depth))
# 
#       for(i in 1:length(x$subtree))
#       {
#         print.Rcknnpart(x$subtree[[i]],depth=parms$depth+1,
#                         beta=parms$beta,
#                         featureNames=parms$featureNames,
#                         mode=parms$mode,
#                         y=parms$y)
#       }
#     }
#   }
# }
