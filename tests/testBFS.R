library(tidyverse)
library(dbrpart)

rm(list=ls())

if(TRUE){

df <- iris

do.call("c",do.call("c",
                    lapply(1:2,function(x){apply(combn(colnames(df)[1:4],x),2,list)})
)) -> combinations

# combinations <- list("Sepal.Length","Petal.Length")

lapply(combinations,function(x){
  if(length(x)>1)
  {
    out <- apply(df[,x],1,paste,collapse=", ")
  }
  else
  {
    out <- df[,x]
  }
  paste("(",out,")",sep="")
}) %>% do.call("cbind",.) -> descriptors

lapply(combinations,function(x){as.matrix(dist(df[,x]))}) -> distanceList
names(distanceList) <- sapply(combinations,paste,collapse=":")


# Code response for Rcpp
target_discrete <- factor(df$Species)
target_continuous <- 0
target_distance <- matrix(0)
mode <- 0

# Swap to zero indexing
target_discrete <- as.numeric(target_discrete)-1
nObs <- length(target_discrete)

# Code metadata about response for Rcpp
classNames <- levels(factor(df$Species))
nCategories <- length(classNames)
classFreq <- as.vector(table(factor(df$Species)))

priors<- classFreq/sum(classFreq)


loss <- 1-diag(nCategories)
acceptableError <- sapply(distanceList,function(x){min(x[x>0],na.rm = TRUE)})/2


controls <- dbrpart.control()

foldList <- rep((1:10)-1, length.out=length(df$Species))
foldList <- foldList[order(runif(length(foldList)))]

} # End bunch of loading and processing


# Debugging the hell out of bubble blower
# RcppParallel::setThreadOptions(numThreads = 1)
# RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())

# sink(file="log.txt",append = FALSE)
fit_dbrpart(in_target_discrete = target_discrete,
            in_target_continuous = target_continuous,
            in_target_distance = target_distance,
            in_factorList = list(),
            in_numericList = list(),
            in_distanceList = distanceList,
            in_acceptable_error=acceptableError,
            in_mode = 0,
            in_stopMethod = controls$stopMethod,
            in_nSims = controls$nSims,
            in_signifLevel = controls$signifLevel,
            in_fold = foldList,
            in_nCategories = nCategories,
            in_classFreq =classFreq,
            in_priors = priors,
            in_loss =loss,
            in_minSplit= controls$minSplit,
            in_minBucket= controls$minBucket,
            in_maxDepth= controls$maxDepth,
            in_alteredPriors= FALSE,
            in_trace= 4,
            in_classMethod= 2, #Classify according to modal value
            in_pCenters= controls$pCenters
            ) -> x
x
# x
# sink()


# x

