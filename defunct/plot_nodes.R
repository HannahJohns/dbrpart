#' Prepare plots showing splits for each node
#'
#'
#'
plot_nodes <- function(tree,target,plotList,featureNames=NULL,featureIDs=NULL,alpha=NULL)
{
  if("Rcknnpart" %in% class(tree))
  {
    alpha <- tree$bestBeta

    targetPlot <- target(tree$featureIDs[tree$tree$groupMembers+1],NULL)
    if(!tree$tree$terminal) # If we split here
    {
      if(alpha<=tree$tree$criticalAlpha) # If we haven't pruned
      {
        # First level plot will only have the target variable distribution
        targetPlot <- target(tree$featureIDs[tree$tree$inTree$groupMembers+1],
                             tree$featureIDs[tree$tree$outTree$groupMembers+1]
        )
      }
    }

    # Construct the values to be contained in this node
    nodePlots <- list(label=tree$tree$label,
                      terminal = tree$tree$terminal,
                      target=targetPlot)


    if(!tree$tree$terminal) # If we split here
    {
      if(alpha<=tree$tree$criticalAlpha) # If we haven't pruned
      {

        inSplit <- list(plotList[[tree$featureNames[tree$tree$bestD+1]]](tree$featureIDs[tree$tree$inTree$groupMembers+1]))
        outSplit <- list(plotList[[tree$featureNames[tree$tree$bestD+1]]](tree$featureIDs[tree$tree$outTree$groupMembers+1]))

        nodePlots <- c(
          nodePlots,
          inPlot = inSplit,
          outPlot = outSplit
        )

      }
    }

    class(nodePlots) <- c("nodePlot",class(nodePlots))

    outList <- list(nodePlots)
    names(outList) <- tree$tree$label

    if(!tree$tree$terminal) # If we split here
    {
      if(alpha<=tree$tree$criticalAlpha) # If we haven't pruned
      {
        outList <- c(
          outList,
          plot_nodes(tree$tree$inTree,target,plotList,featureNames=tree$featureNames,featureIDs=tree$featureIDs,alpha=alpha),
          plot_nodes(tree$tree$outTree,target,plotList,featureNames=tree$featureNames,featureIDs=tree$featureIDs,alpha=alpha)
        )
      }
      else
      {
      }
    }
    else
    {
    }

  }
  else
  {

    targetPlot <- target(tree$featureIDs[tree$tree$groupMembers+1],NULL)

    if(!tree$terminal) # If we split here
    {
      if(alpha<=tree$criticalAlpha) # If we haven't pruned
      {
        targetPlot <- target(featureIDs[tree$inTree$groupMembers+1],
                             featureIDs[tree$outTree$groupMembers+1]
        )
      }
    }

    # Construct the values to be contained in this node
    nodePlots <- list(label=tree$label,
                      terminal = tree$terminal,
                      target=targetPlot)


    if(!tree$terminal) # If we split here
    {
      if(alpha<=tree$criticalAlpha) # If we haven't pruned
      {

        inSplit <- list(plotList[[featureNames[tree$bestD+1]]](featureIDs[tree$inTree$groupMembers+1]))
        outSplit <- list(plotList[[featureNames[tree$bestD+1]]](featureIDs[tree$outTree$groupMembers+1]))

        nodePlots <- c(
          nodePlots,
          inPlot = inSplit,
          outPlot = outSplit
        )

      }
    }

    class(nodePlots) <- c("nodePlot",class(nodePlots))

    outList <- list(nodePlots)
    names(outList) <- tree$label

    if(!tree$terminal) # If we split here
    {
      if(alpha<=tree$criticalAlpha) # If we haven't pruned
      {
        outList <- c(
          outList,
          plot_nodes(tree$inTree,target,plotList,featureNames=featureNames,featureIDs=featureIDs,alpha=alpha),
          plot_nodes(tree$outTree,target,plotList,featureNames=featureNames,featureIDs=featureIDs,alpha=alpha)
        )
      }
    }

  }
  return(outList)
}

plot.nodePlot <- function(x,...){
  if(x$terminal)
  {
    invisible(plot(x$target))
  }
  else
  {
    invisible(plot(cowplot::plot_grid(x$inPlot,
                                      cowplot::plot_to_gtable(
                                        x$target +
                                        ggplot2::labs(title=sprintf("Node %s",x$label)) +
                                        ggplot2::theme(title=element_text(hjust = 0.5))
                                      ),
                                      x$outPlot,
                                      nrow=1,
                                      ...)
                   ))
  }
}
