#' Cluster record IDs into terminal nodes
#'
#' Produce clusters based on terminal nodes
#' of tree
#' 
#' @param x
#' 
#' @returns A matrix containing the record number and the terminal node
#'          it corresponds to
#'
#' @export
cluster <- function(x,step=0,...)
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
  
  if(x$tree$terminal[parms$label] | x$tree$pruneAtStep[parms$label] <= (parms$pruneStep))
  {
    out <- cbind(x$tree$groupMembers[[parms$label]]+1,parms$label)
  }
  else
  {
    out <- rbind(
              cluster(x,label=x$tree$inChild[parms$label]+1,pruneShift = parms$pruneStep),
              cluster(x,label=x$tree$outChild[parms$label]+1, pruneShift = parms$pruneStep)
    )
  }
  
  return(out)
  
  # 
  # parms <- list(...)
  # # Recursively traverse the tree and group by terminal nodes
  # if("Dbrpart" %in% class(x))
  # {
  #   if(is.null(parms$beta)){parms$beta <- x$bestBeta}
  # 
  #   result <- cluster(x$tree,beta=parms$beta)
  # 
  #   # Convert from position IDs to row names
  #   # and order the results appropriately
  # 
  #   out <- result[,2]
  #   names(out) <- x$featureIDs[result[,1]]
  # 
  #   out <- out[match(x$featureIDs, names(out))]
  # 
  #   return(out)
  # 
  # } else {
  # 
  #   if(x$terminal | (parms$pruneAtStep >= step))
  #   {
  #     
  #   }
  #   else
  #   {
  #     out <- rbind(cluster(x$subtree[[1]],beta=parms$beta),
  #                  cluster(x$subtree[[2]],beta=parms$beta))
  #   }
  #   return(out)
  # }
}

