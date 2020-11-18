treeGrowth <- function(x,depth=0)
{

  if(sum(x$a <= 0)>1)
  {
    x$nTerm <- x$nTerm[-1]
  }
  
  df <- data.frame(beta=rev(x$b),
                   nTerminal=rev(x$nTerm),
                   risk=rev(x$risk),
                   cv=rev(x$cv_risk))

  df <- unique(df)

  # If we have two competing trees of the same size,
  # only include the one with smallest risk
  df <- do.call("rbind",by(df,df$nTerminal,function(y){

    if(nrow(y)>1)
    {
      return(y[which(y$risk==min(y$risk)),])
    }
    else
    {
      return(y)
    }
  }))

  # If we have two competing trees of the same risk,
  # only include the one with smallest size
  df <- do.call("rbind",by(df,df$risk,function(y){

    if(nrow(y)>1)
    {
      return(y[which(y$nTerminal==min(y$nTerminal)),])
    }
    else
    {
      return(y)
    }
  }))

  df <- df[order(df$nTerminal),]

  df$beta[1] <- Inf
  df$relative.risk <- df$risk/df$risk[1]
  df$relative.risk_cv <- df$cv/df$cv[1]


  beta <- df$beta[which(df$cv==min(df$cv))+depth]
  nTerminal <- df$nTerminal[which(df$cv==min(df$cv))+depth]
  relative.risk <- df$relative.risk[which(df$cv==min(df$cv))+depth]
  relative.risk_cv <- df$relative.risk_cv[which(df$cv==min(df$cv))+depth]
  
  tieBreak <- which(beta==max(beta))
  beta <- beta[tieBreak]
  nTerminal <- nTerminal[tieBreak]
  relative.risk <- relative.risk[tieBreak]
  relative.risk_cv <- relative.risk_cv[tieBreak]
  
  
  out <- list(frame=df,
              beta=beta,
              nTerminal=nTerminal,
              relative.risk=relative.risk,
              relative.risk_cv=relative.risk_cv)

  class(out) <- c("TreeGrowth",class(out))

  return(out)
}

plot.TreeGrowth <- function(x,...)
{
  parms <- list(...)

  if(is.null(parms$direction))
  {
    parms$direction <- "risk"
  }

  plot(1:nrow(x$frame),x$frame$relative.risk,axes=FALSE,xlab="",ylab="",col="blue",type="b",
       ylim=c(min(x$frame$relative.risk,x$frame$relative.risk_cv),1))
  lines(1:nrow(x$frame),x$frame$relative.risk_cv,col="red")
  points(1:nrow(x$frame),x$frame$relative.risk_cv,col="red")
  
  if(parms$direction =="risk")
  {
    axis(2,x$frame$relative.risk,round(x$frame$relative.risk,2),las=1)
  }

  if(parms$direction =="rsquare")
  {
    axis(2,x$frame$relative.risk,round(1-x$frame$relative.risk,2),las=1)
  }

  axis(1,1:nrow(x$frame),round(x$frame$beta,3),las=1,line=1)
  mtext("Beta",1,line=1,at=0.8)
  axis(1,1:nrow(x$frame),x$frame$nTerminal,las=0,line=3)
  mtext("nTerm",1,line=3,at=0.8)
  abline(v=which(x$frame$beta==x$beta),col="red",lty=3)
  

  legend("topright",c("Model","Model under CV"),lty=1,col=c("blue","red"))

  if(parms$direction =="risk")
  {
    title("Relative Risk of Tree")
  }
  else if(parms$direction =="rsquare")
  {
    title("r-square of Tree")
  }
}

# Plot and first axis:
plot(1:10,1:10,bty="n",col="red",pch=16,axes=FALSE,xlab="",ylab="")
axis(2,0:11,las=1)
axis(1,0:11,line=1,col="red",col.ticks="red",col.axis="red")

# Secondary points and axis:
points(rnorm(10,50,20)/10, rnorm(10,5,2),pch=16, col="blue" )
axis(1,0:11,labels=0:11*10,line=3,col="blue",col.ticks="blue",col.axis="blue")
mtext("Label 2",1,line=3,at=0.2,col="blue")
