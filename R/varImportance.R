varImportance <- function(x){

  recurseImportance <- function(y){
    out <- c(var=y$bestD,
             reduction=y$probability*y$impurity -
                       (
                         y$inTree$probability*y$inTree$impurity+
                         y$outTree$probability*y$outTree$impurity
                       )
             )

    # Only consider splits
    if(!y$inTree$terminal)
    {
      out <- rbind(out,recurseImportance(y$inTree))
    }

    if(!y$outTree$terminal)
    {
      out <- rbind(out,recurseImportance(y$outTree))
    }

    rownames(out) <- NULL
    return(out)
  }

  nodeImportance <- recurseImportance(x$tree)

  out <- tapply(nodeImportance[,"reduction"],x$featureNames[nodeImportance[,"var"]+1],sum)
  return(out/sum(out))
}
