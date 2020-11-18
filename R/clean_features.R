clean_features <- function(x,DAG=NULL,epsilon=0.0001,combine=NULL)
{

  changeList <- list(inflation=rep(0,length(x)),
                     scale=rep(0,length(x)),
                     keep={},
                     order=NA)

  # Name each feature
  if(is.null(names(x)))
  {
    names(x) <- as.character(1:length(x))
  }

  #Get record label vectors for all features
  nameList <- lapply(x,function(y){
    if(class(y)=="dist")
    {
      return(attr(y,"Labels"))
    }
    else if(class(y)=="matrix")
    {
      if(is.null(rownames(y)) | is.null(colnames(y)))
      {
        return(NULL)
      }
      if(prod(rownames(y)==rownames(y))==1)
      {
        return(rownames(y))
      }
      else
      {
        stop("Row and Column names don't match")
      }
    }
  })

  # Cast everything to matrices
  x <- lapply(x,function(y){as.matrix(y)})

  # Check that name vectors are consistent
  if(sum(sapply(nameList,is.null))==length(nameList))
  {
    warning("Recommend providing labels or row names for all records")
  }
  else
  {
    if(length(unique(nameList))>1)
    {
      if(length(unique(lapply(nameList,sort)))==1)
      {
        warning("Feature lists do not have consistent ordering for labels. Re-arranging records for all features")
        changeList$order <- rownames(x[[1]])
        x<- lapply(x,function(y){
              y[rownames(y)[order(match(rownames(y),changeList$order))],
                rownames(y)[order(match(rownames(y),changeList$order))]
                ]
            })
      }
      else
      {
        stop("labels are not consistent across features")
      }
    }
  }


  nrows <- unique(sapply(x,function(y){nrow(y)}))
  if(length(nrows)>1)
  {
    stop("Matrices are of different sizes")
  }

  if(prod(sapply(x,function(y){nrow(y)==ncol(y)}))==0)
  {
    stop("Matrices are not all rectangular")
  }

  # Off-diagonals must be stricly positive. Inflate distances if needed
  for(i in 1:length(x))
  {
    if(min(x[[i]][row(x[[i]])!=col(x[[i]])],na.rm = TRUE)==0)
    {
      x[[i]] <- x[[i]]+epsilon
      changeList$inflation[i] <- epsilon
    }
  }

  # Rescale distance matrices to have offdiagonal average of 1
  for(i in 1:length(x))
  {
    changeList$scale[i] <- mean(x[[i]][row(x[[i]])!=col(x[[i]])],na.rm=TRUE)
    x[[i]] <- x[[i]]/changeList$scale[i]
  }

  # Get list of rows/columns which are not missing in every matrix
  # Could potentially do some fancy things in here re: correcting missing values
  # but that's a much harder problem
  keepList <- 1:nrows
  for(i in 1:length(x))
  {
   # Allow for at least one nonzero entry (diagonal)
    keepList <- base::intersect(keepList,which(!apply(x[[i]],1,function(y){sum(is.na(y))>=(nrow(x[[i]])-1)})))
  }
  changeList$keep <- keepList

  # If no DAG is supplied, assume no link between predictors
  if(is.null(DAG))
  {
    DAG <- diag(nrow=length(x))
  }

  # Get penalty adjustments

  out <- list(featureList = x,
              DAG = DAG,
              changeList = changeList)

  class(out) <- c("cleaned.distanceList",class(out))

  return(out)

}

summary.cleaned.distanceList <- function(x,...)
{
  cat(sprintf("Cleaned %d distance matrices, resulting in:\n",length(x$featureList)))

  if(sum(x$changeList$inflation>0)>0)
  {
    cat(sprintf("\t Distance inflated on features %s\n",
                paste(names(x$featureList)[x$changeList$inflation>0],
                      collapse=", ")))
  }
  else
  {
    cat("\t No distance inflation applied\n")
  }

  cat("\t Features scale reduced by factor of\n")
  cat(sprintf("\t\t %s %f\n",names(x$featureList),x$changeList$scale))

  nDropped <- nrow(x$featureList[[1]]) - length(x$changeList$keep)
  cat(sprintf("\t %d records identified with missing features\n",nDropped))
  cat("\n")
}


export.cleaned.distanceList <- function(x,...)
{
  lapply(x$featureList,function(D){
    as.dist(D[x$changeList$keep,x$changeList$keep])
  })
}


