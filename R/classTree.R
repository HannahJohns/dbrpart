classTree <- function(x)
{
  classif <- x$class
  colnames(classif) <- paste("L",1:ncol(classif),sep="")

  # Strip excess classification depths
  i=1
  while(!prod(as.character(classif[,i])==as.character(classif[,i+1])))
  {
    i=i+1
  }
  classif <- as.data.frame(classif[,1:i])

  classif$y <- x$y
  if(!is.null(x$classNames))
  {
    classif$y <- x$classNames[classif$y+1]
  }

  if(!is.null(x$featureIDs))
  {
    classif$id=x$featureIDs
  }

  class(classif) <- c("Rcknpart.ClassTree",class(classif))

  return(classif)
}

plot.Rcknpart.ClassTree <- function(x,...)
{
  if(!requireNamespace("treemap",quietly = TRUE))
  {
   stop("treemap package needed for this function to work. Please install it.")
  }

  parms <- list(...)
  if(is.null(parms$depth))
  {
    depth <- length(grep("L",colnames(x)))
  }

  if(is.null(parms$palette))
  {
    parms$palette <- rep("Blues",length(levels(x$y)))
  }

  # Need to summarise these classifications before we can use it

  classifSummary <- {}

  for(class in unique(x$y))
  {
    classifSummary <- cbind(classifSummary,
                            tapply(x$y,
                                   x[,paste("L",depth,sep="")],
                                   function(y){sum(y==class)}))
  }
  colnames(classifSummary) <- unique(x$y)
  classifSummary <- as.data.frame(classifSummary)

  count <- apply(classifSummary,1,sum)
  for(class in unique(x$y))
  {
    classifSummary[,class] <- classifSummary[,class]/count
  }
  classifSummary <- cbind(classifSummary,count)

  # If in the odd case COUNT is a level, we can't use it.
  # Find some other variable name we can use instead.
  if("COUNT" %in% unique(x$y))
  {
    i <- 1
    while(sprintf("COUNT%d",i) %in% unique(x$y))
    {
      i <- i+1
    }
    countName <- sprintf("COUNT%d",i)
  } else {
    countName <- "COUNT"
  }
  colnames(classifSummary)[ncol(classifSummary)] <- countName

  classifSummary <- cbind(classifSummary,rownames(classifSummary))
  colnames(classifSummary)[ncol(classifSummary)] <- paste("L",depth,sep="")

  classifSummary <- base::merge(classifSummary,unique(x[,grep("L",colnames(x))]))
  classifSummary <- classifSummary[,c(paste("L",1:depth,sep=""),countName,unique(x$y))]

  out <- list()

  i <- 0
  for(class in unique(x$y))
  {
    i <- i+1

    out <- c(out,
             treemap::treemap(classifSummary,
                              index=paste("L",1:depth,sep=""),
                              vSize = countName,vColor = class,type="value",
                              n=10, palette = parms$palette[i],
                              fontsize.labels=c(rep(0,depth-1),10),title = "")
    )
  }

  return(invisible(out))
}
