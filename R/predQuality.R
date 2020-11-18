

predQuality <- function(fit,index){
  
  qualityTable <- {}
  
  for(maxmRS in 0:5){
    
    
    auc <- as.numeric(roc(
      as.numeric(fit$y)-1 <= maxmRS,
      apply(fit$details[[index]],1,function(x){sum(x[1:(maxmRS+1)])})
    )$auc)
    
    confusion <- table(
      as.numeric(fit$y)-1 <= maxmRS,
      fit$prediction[,index]<=maxmRS
    )
    
    if(ncol(confusion)==1 & colnames(confusion)[1]=="FALSE")
    {
      confusion <- cbind(confusion,0)
    }
    
    if(ncol(confusion)==1 & colnames(confusion)[1]=="TRUE")
    {
      confusion <- cbind(0,confusion)
    }
    
    
    
    quality <-c(
      auc=auc,
      spec=confusion[1,1]/sum(confusion[1,]), sens=confusion[2,2]/sum(confusion[2,]),
      npv=confusion[1,1]/sum(confusion[,1]), ppv=confusion[2,2]/sum(confusion[,2])
    )
    
    cv_auc <- as.numeric(roc(
      as.numeric(fit$y)-1 <= maxmRS,
      apply(fit$cv_details[[index]],1,function(x){sum(x[1:(maxmRS+1)])})
    )$auc)
    
    cv_confusion <- table(
      as.numeric(fit$y)-1 <= maxmRS,
      fit$cv_prediction[,index]<=maxmRS
    )
    
    if(ncol(cv_confusion)==1 & colnames(cv_confusion)[1]=="FALSE")
    {
      cv_confusion <- cbind(cv_confusion,0)
    }
    
    if(ncol(cv_confusion)==1 & colnames(cv_confusion)[1]=="TRUE")
    {
      cv_confusion <- cbind(0,cv_confusion)
    }
    
    
    cv_quality=c(
      cv_auc=cv_auc,
      cv_spec=cv_confusion[1,1]/sum(cv_confusion[1,]), cv_sens=cv_confusion[2,2]/sum(cv_confusion[2,]),
      cv_npv=cv_confusion[1,1]/sum(cv_confusion[,1]), cv_ppv=cv_confusion[2,2]/sum(cv_confusion[,2])
    )
    
    qualityTable <- rbind(qualityTable,c(maxmRS = maxmRS,quality,cv_quality))
  }
  
  qualityTable <- qualityTable[,c("maxmRS",
                                  "auc","cv_auc",
                                  "spec","cv_spec",
                                  "sens","cv_sens",
                                  "npv","cv_npv",
                                  "ppv","cv_ppv")]
  
  qualityTable <- as.data.frame(qualityTable)
  
  return(qualityTable)
}
