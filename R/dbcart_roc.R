dbcart_roc <- function(x,beta=NULL)
{

  out <- data.frame()

  if(is.null(beta)){beta <- x$bestBeta}

  classProbs <- predict(x,beta=beta)
  classProbs_cv <- x$classProbs_cv[[max(which(x$betas<=beta))]]

  classProbs <- as.matrix(classProbs[,-ncol(classProbs)])
  classProbs_cv <- classProbs_cv[,-ncol(classProbs)]

  for(i in 1:(length(x$classNames)-1)){

    cut <- sprintf("( %s ) | ( %s )",
                   paste(x$classNames[1:i], collapse="-"),
                   paste(x$classNames[(i+1):length(x$classNames)],collapse="-"))

    broc <- pROC::roc(x$y<=(i-1),
                apply(classProbs,1,function(x){sum(x[1:i])}))

    broc$auc

    proc <- pROC::roc(x$y<=(i-1),
                apply(classProbs_cv,1,function(x){sum(x[1:i])}))

    proc$auc

    threshold <- with(broc,thresholds[sensitivities+specificities == max(sensitivities+specificities)])

    rbind(
      diag(prop.table(
        table(
          paste("pred",apply(classProbs,1,function(x){sum(x[1:i])})>threshold),
          paste("real",x$y<=(i-1))
        ),
        margin=2
      )),
      diag(
        prop.table(
          table(
            paste("pred",apply(classProbs_cv,1,function(x){sum(x[1:i])})>threshold),
            paste("real",x$y<=(i-1))
          ),
          margin=2
        )
      )) -> sensspec
    colnames(sensspec) <- c("Specificity","Sensitivity")
    rownames(sensspec) <- c("Model","CV")

    rbind(
      diag(prop.table(
        table(
          paste("pred",apply(classProbs,1,function(x){sum(x[1:i])})>threshold),
          paste("real",x$y<=(i-1))
        ),
        margin=1
      )),
      diag(
        prop.table(
          table(
            paste("pred",apply(classProbs_cv,1,function(x){sum(x[1:i])})>threshold),
            paste("real",x$y<=(i-1))
          ),
          margin=1
        )
      )) -> ppvnpv
    colnames(ppvnpv) <- c("npv","ppv")
    rownames(ppvnpv) <- c("Model","CV")


    out <- rbind(out,
                 data.frame(
                   cut = cut,
                   roc.model = as.numeric(broc$auc),
                   roc.cv = as.numeric(proc$auc),
                   sens.model = sensspec["Model","Sensitivity"],
                   sens.cv = sensspec["CV","Sensitivity"],
                   spec.model = sensspec["Model","Specificity"],
                   spec.cv = sensspec["CV","Specificity"],
                   ppv.model = ppvnpv["Model","ppv"],
                   ppv.cv = ppvnpv["CV","ppv"],
                   npv.model = ppvnpv["Model","npv"],
                   npv.cv = ppvnpv["CV","npv"]
                 )
                 )


  }

  class(out) <- c("Dbcartroc",class(out))

  return(out)
}

print.Dbcartroc <- function(x,...)
{
  for(i in 1:nrow(x))
  {
    cat(as.character(x[i,"cut"]),"\n")

    vals <- matrix(c(x[i,"roc.model"],x[i,"sens.model"],x[i,"spec.model"],x[i,"ppv.model"],x[i,"npv.model"],
                     x[i,"roc.cv"],   x[i,"sens.cv"],   x[i,"spec.cv"],   x[i,"ppv.cv"],   x[i,"npv.cv"]),
                   nrow=2,byrow = TRUE)

    colnames(vals) <- c("roc","sensitivity","specificity","ppv","npv")
    rownames(vals) <- c("Model","CV")

    print(vals)
  }

}


