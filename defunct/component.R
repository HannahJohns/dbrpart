component <- function(x,trim_k=0.6)
{
  features <- x$input$features
  outcome <- x$input$outcome
  ks <- x$input$ks
  metaData <- x$input$metaData
  nFolds <- x$input$nFolds
  evalMethod <- x$input$evalMethod
  somaOptions <- x$input$somaOptions

  if(length(ks>1))
  {
    ks <- x$wcknn$kList[x$wcknn$kList[,2]>quantile(x$wcknn$kList[,2],0.5),1]
  }

  featureNames <- names(features)
  if(is.null(featureNames))
  {
    featureNames <- as.numeric(1:length(features))
  }

  somaOptions$verbose <- somaOptions$verbose-1
  altModels <- list()
  # Get the effect of removing each feature on model quality and compare outcomes
  for(i in 1:length(features))
  {
    if(somaOptions$verbose>=0)
    {
      cat(sprintf("\nTesting Feature %s\n",featureNames[i]))
    }

     X <- cbind(0,
               matrix(rnorm(somaOptions$popSize*((length(features)-2))),
                      ncol=(length(features)-2)
               )
    )

    fit <- grandDistance_fit_k(D = features[-i], X = X,outcome = outcome,
                               metaData = metaData,
                               k = ks,nFolds = nFolds,evalMethod = evalMethod,
                               somaOptions = somaOptions)
    altModels[[i]] <- fit
  }

  dValue <- sapply(altModels,function(y,eval){eval-y$bestEval},eval=x$wcknn$bestEval)
  names(dValue) = featureNames

  out <- list(wcknn=x,dValue=dValue,altModels=altModels)

  class(out) <- c(class(out),"component")

  return(out)
}

plot.component <- function(x)
{
  # Compare training results for each submodel

  kData <- c(range(x$wcknn$wcknn$kList[,1],finite=TRUE),
             range(x$wcknn$wcknn$kList[,2],finite=TRUE))
  yMin <- min(kData[3])
  yMax <- max(kData[4])
  kData <- do.call("rbind",
                   lapply(x$altModels,
                          function(x){c(range(x$kList[,1],finite=TRUE),
                                        range(x$kList[,2],finite=TRUE))}
                   ))
  kMin <- min(kData[,1])
  kMax <- max(kData[,2])
  yMin <- min(yMin,kData[,3])
  yMax <- max(yMax,kData[,4])

  plot.new()
  for(i in 0:length(interpret$altModels))
  {
    par(new=TRUE)
    if(i==0)
    {
      plot(x$wcknn$wcknn$kList[,1],
           x$wcknn$wcknn$kList[,2],
           xlim = c(kMin,kMax),
           ylim = c(yMin,yMax),
           pch=2,
           col=i+1,
           xlab="k",
           ylab=x$wcknn$input$evalMethod$method)
    }
    else
    {
      plot(x$altModels[[i]]$kList[,1],
           x$altModels[[i]]$kList[,2],
           type="l",
           xlim = c(kMin,kMax),
           ylim = c(yMin,yMax),
           pch=3,
           col=i+1,axes = FALSE,sub = "",xlab = "",ylab = "")
    }
  }
  legend("bottomleft",
         lty=c(NA,rep(1,length(interpret$altModels))),
         pch=c(2,rep(NA,length(interpret$altModels))),
         col=(0:length(interpret$altModels))+1,
         legend=c("Full Model",paste("Drop",names(x$wcknn$input$features))))

  # Compare predictions

}



