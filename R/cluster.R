#' Cluster record IDs into terminal nodes
#'
#' Produce clusters based on terminal nodes
#' of tree
#'
#' @export
cluster <- function(x,...)
{
  parms <- list(...)
  # Recursively traverse the tree and group by terminal nodes
  if("Rcknnpart" %in% class(x))
  {
    if(is.null(parms$beta)){parms$beta <- x$bestBeta}

    result <- cluster(x$tree,beta=parms$beta)

    # Convert from position IDs to row names
    # and order the results appropriately

    out <- result[,2]
    names(out) <- x$featureIDs[result[,1]]

    out <- out[match(x$featureIDs, names(out))]

    return(out)

  } else {

    if(x$terminal | (parms$beta >= x$criticalAlpha))
    {
      out <- cbind(x$groupMembers+1,x$label)
    }
    else
    {
      out <- rbind(cluster(x$subtree[[1]],beta=parms$beta),
                   cluster(x$subtree[[2]],beta=parms$beta))
    }
    return(out)
  }
}

