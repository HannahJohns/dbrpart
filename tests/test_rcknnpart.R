rm(list=ls())

library(dbrpart)

# Datasets, see https://machinelearningmastery.com/machine-learning-datasets-in-r/
library(mlbench)

set.seed(NULL)
library(tidyverse)
library(geosphere)
library(rpart)
#library(hannahTools)





##########################################################################
#
#                         DISTANCE-BASED REGRESSION TREES
#
#        Note - objective matrix will be dist(dbind(zn, tax, cmedv))
#
#


rm(list=ls())

data(BostonHousing2)

  

  # straightLineDistance <- matrix(0,nrow=nrow(BostonHousing2),ncol=nrow(BostonHousing2))
  # for(i in 1:(nrow(BostonHousing2)-1))
  # {
  #   print(i)
  #   for(j in (i+1):nrow(BostonHousing2))
  #   {
  #     d <- distVincentyEllipsoid(BostonHousing2[i,c("lat","lon")],BostonHousing2[j,c("lat","lon")])
  #     straightLineDistance[i,j] <- d
  #     straightLineDistance[j,i] <- d
  #   }
  # }
  # 
  
  featureList <- with(BostonHousing2,
                      list(
                        lat=dist(lat),
                        lon=dist(lon)
                        # spatial = as.dist(straightLineDistance)
                      )
  )
  featureList <- lapply(featureList,as.dist)
  
  value <- BostonHousing2$cmedv
  
  range(BostonHousing2$lat)
  range(BostonHousing2$lon)
  
  
  value <- dist(BostonHousing2$cmedv)

outDB <- dbrpart(target=value,
                 distanceList = featureList,
                 controls = dbrpart.control(
                   trace=999,
                   minSplit = 160,
                   minBucket = 80,
                   pCenters = 0,
                   stopMethod = 1
                 ))





plot(treeGrowth(out))

out$betas
out$risks

treeGrowth(out)


out$tree$risk
sum(sapply(out$tree$subtree,function(x){x$probability*x$risk}))

quasiregression_start <- Sys.time()

quasiregression_stop <- Sys.time()
quasiregression_stop-quasiregression_start

plot(treeGrowth(out))

















##################################################

rm(list=ls())

df <- iris

lapply(1:4,function(x){dist(df[,x])}) -> D_list
names(D_list) <- colnames(df)[1:4]


do.call("c",do.call("c",
                    lapply(1:2,function(x){apply(combn(0:3,x),2,list)})
)) -> combinations


# First off just test that printing is fine
# 
# out <- dbrpart(target=df$Species,distanceList = D_list, distanceCombinations = combinations,
               # controls = dbrpart.control(classMethod = "mode",pCenters = 0,trace = 3))


sink("log.txt")
out <- dbrpart(target=df$Species,distanceList = D_list,
                controls = dbrpart.control(classMethod = "mode",pCenters = 0,trace = 3,stopMethod = 1))
sink()

print(out)

treeGrowth(out)


plot(treeGrowth(out))

table(out$y,out$prediction[,1])

#####################################################

# Parallelisation tests


RcppParallel::setThreadOptions(numThreads = 1)
sink(file="log.txt")
out <- dbrpart(target=df$Species,distanceList = D_list,
               controls = dbrpart.control(classMethod = "mode",pCenters = 0,trace = 999,minBucket = 20,minSplit = 60))
sink()
table(out$y,out$prediction[,1])

RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())
sink(file="log2.txt")
out2 <- dbrpart(target=df$Species,distanceList = D_list,
                controls = dbrpart.control(classMethod = "mode",pCenters = 0,trace = 999,minBucket = 20,minSplit = 60))
sink()
table(out2$y,out2$prediction[,1])


plot(treeGrowth(out2))



rbind(
cbind(out$tree$groupMembers[[3]]+1,out$tree$label[3]),
cbind(out$tree$groupMembers[[4]]+1,out$tree$label[4]),
cbind(out$tree$groupMembers[[5]]+1,out$tree$label[5])
) -> termNodesReal

table(termNodesReal[order(termNodesReal[,1]),2],out$classTree[,1])

out$tree$groupMembers[[4]]+1
out$tree$groupMembers[[5]]+1
#####################################

# Debug tests

rm(list=ls())


# "Pen and Paper" calculations

df <- iris

lapply(1:4,function(x){as.matrix(dist(df[,x]))}) -> D_list
y <- as.numeric(df$Species)
priors <- as.numeric(prop.table(table(y)))
members <- 1:length(y)


do.call("rbind",lapply(1:4,function(feature){
do.call("rbind",lapply(members,function(c){
  radius <- sort(unique(D_list[[feature]][c,members]))
  utility <- sapply(radius,function(r){
    inList <- members[D_list[[feature]][c,members] <= r]
    outList <- members[D_list[[feature]][c,members] > r]

    inTab <- sapply(1:3,function(i){sum(y[inList]==i)})
    outTab <- sapply(1:3,function(i){sum(y[outList]==i)})

    if(sum(inTab) < 7 | sum(outTab) < 7 )
    {
      return(Inf)
    }
    
    inP <- sum(priors*inTab/50)
    outP <- sum(priors*outTab/50)  
            
    inPclass <- (priors*inTab/50) / inP
    outPclass <- (priors*outTab/50) / outP
    
    gini_in <- sum(sapply(1:3,function(i){sapply((1:3)[-i],function(j){inPclass[i] * inPclass[j]})}))
    gini_out <- sum(sapply(1:3,function(i){sapply((1:3)[-i],function(j){outPclass[i] * outPclass[j]})}))
    
    impurity <- inP * gini_in + outP * gini_out 
    
    return(impurity)
  })
  
  return(c(c,max(radius[utility==min(utility)])[1],min(utility)[1]))
})) -> bubbles

bestBubbles <- bubbles[bubbles[,3]==min(bubbles[,3]),]
bestBubbles <- bestBubbles[bestBubbles[,2]==max(bestBubbles[,2]),]
if("matrix" %in% class(bestBubbles))
{
  bestBubbles <- bestBubbles[1,]
}
return(c(feature,bestBubbles))
})) -> splits



inGroup <- members[D_list[[splits[3,1]]][splits[3,2],members] <= splits[3,3]]
members <- inGroup

do.call("rbind",lapply(1:4,function(feature){
  do.call("rbind",lapply(members,function(c){
    radius <- sort(unique(D_list[[feature]][c,members]))
    utility <- sapply(radius,function(r){
      inList <- members[D_list[[feature]][c,members] <= r]
      outList <- members[D_list[[feature]][c,members] > r]
      
      inTab <- sapply(1:3,function(i){sum(y[inList]==i)})
      outTab <- sapply(1:3,function(i){sum(y[outList]==i)})
      
      if(sum(inTab) < 7 | sum(outTab) < 7 )
      {
        return(Inf)
      }
      
      inP <- sum(priors*inTab/50)
      outP <- sum(priors*outTab/50)  
      
      inPclass <- (priors*inTab/50) / inP
      outPclass <- (priors*outTab/50) / outP
      
      gini_in <- sum(sapply(1:3,function(i){sapply((1:3)[-i],function(j){inPclass[i] * inPclass[j]})}))
      gini_out <- sum(sapply(1:3,function(i){sapply((1:3)[-i],function(j){outPclass[i] * outPclass[j]})}))
      
      impurity <- inP * gini_in + outP * gini_out 
      
      return(impurity)
    })
    
    return(c(c,max(radius[utility==min(utility)])[1],min(utility)[1]))
  })) -> bubbles
  
  bestBubbles <- bubbles[bubbles[,3]==min(bubbles[,3]),]
  bestBubbles <- bestBubbles[bestBubbles[,2]==max(bestBubbles[,2]),]
  if("matrix" %in% class(bestBubbles))
  {
    bestBubbles <- bestBubbles[1,]
  }
  return(c(feature,bestBubbles))
})) -> splits



inGroup <- members[D_list[[splits[3,1]]][splits[3,2],members] <= splits[3,3]]
outGroup <- members[D_list[[splits[3,1]]][splits[3,2],members] > splits[3,3]]










##############################################


RcppParallel::setThreadOptions(numThreads = 1)



df <- iris

lapply(1:4,function(x){dist(df[,x])}) -> D_list
names(D_list) <- colnames(df)[1:4]

out <- dbrpart(target=df$Species,distanceList = D_list,
               controls = dbrpart.control(classMethod = "mode",pCenters = 1,trace = 999,minBucket = 7,minSplit = 30))

table(out$y,out$prediction[,1])

x <- rpart(Species~.,data=df)

table(df$Species,apply(predict(x),1,function(x){(1:3)[x==max(x)]}))

{
sink(file="log.txt")

df <- iris

lapply(1:4,function(x){dist(df[,x])}) -> D_list
names(D_list) <- colnames(df)[1:4]

out <- dbrpart(target=df$Species,distanceList = D_list,
               controls = dbrpart.control(classMethod = "mode",pCenters = 1,trace = 999,minBucket = 7,minSplit = 30))
sink()
}





plot(treeGrowth(out))

out$tree$depths
sapply(out$tree$groupMembers,length)

tapply(
sapply(out$tree$groupMembers,length),
out$tree$depths,
sum)



out$tree$inChild
out$tree$outChild

intersect(out$tree$groupMembers[[4]],out$tree$groupMembers[[5]])

out$tree$bestSplitIndex
out2$tree$bestSplitIndex

out$tree$radius





out$tree$center
out$tree$radius
out$tree$radius

table(out$y,out$prediction[,1])

out$tree$bestSplitIndex
out$tree$center
out$tree$radius
out$tree$inChild
out$tree$outChild

# Root
as.matrix(D_list[[1]])[1,]<0.6

# In child should be pure
out$tree$groupMembers[[2]]+1

# Out child should need a split still
out$tree$groupMembers[[3]]+1

table(out$y[out$tree$groupMembers[[3]]+1],as.matrix(D_list[[4]])[51,out$tree$groupMembers[[3]]+1]<0.35)

plot(treeGrowth(out),direction="risk")


table(apply(predict(rpart::rpart(Species~.,df)),1,function(x){(1:3)[x==max(x)]}),out$prediction[,1])

do.call("c",do.call("c",
                    lapply(1:2,function(x){apply(combn(0:3,x),2,list)})
)) -> combinations

sink("log.txt")
out <- dbrpart(target=df$Species,distanceList = D_list,distanceCombinations = combinations,
               controls = dbrpart.control(classMethod = "mode",pCenters = 1,trace = 5,minBucket = 7,minSplit = 30))
sink()

table(apply(predict(rpart::rpart(Species~.,df)),1,function(x){(1:3)[x==max(x)]}),out$prediction[,1])

{
################################################################################
#
#                         CLASSIFICATION TREES

rm(list=ls())

df <- iris

lapply(1:4,function(x){dist(df[,x])}) -> D_list
names(D_list) <- colnames(df)[1:4]


do.call("c",do.call("c",
                    lapply(1:2,function(x){apply(combn(0:3,x),2,list)})
)) -> combinations

out <- dbrpart(target=df$Species,distanceList = D_list,
               controls = dbrpart.control(classMethod = "mode",pCenters = 0,trace = 5,minBucket = 40,minSplit = 60))


out2 <- dbrpart(target=df$Species,distanceList = D_list,distanceCombinations = combinations,
               controls = dbrpart.control(classMethod = "mode",pCenters = 0,trace = 5,minBucket = 40,minSplit = 60))

plot(treeGrowth(out),direction="risk")


table(out$y,out$prediction[,1])
table(out2$y,out2$prediction[,1])














out <- dbrpart(target=df$Species,distanceList = D_list,distanceCombinations = as.list(0:3),
               controls = dbrpart.control(classMethod = "mode",pCenters = 0,trace = 4,minBucket = 7,minSplit = 30))




out <- dbrpart(target=df$Species,distanceList = D_list,
               controls = dbrpart.control(classMethod = "mode",pCenters = 1,trace = 999,minBucket = 20,minSplit = 40))
out$tree$pruneAtStep

x <- out$classTree
i <- ncol(x)
while(i>1){
  print(table(x[,i],x[,i-1]))
  i <- i-1
}

out

treeGrowth(out)

plot(treeGrowth(out),direction="risk")
plot(treeGrowth(out),direction="rsquare")

print(out)
cluster(out)

table(df$Species,predict(out)$classif)

# Permutation stoppage is broken
out <- rcknnpart(target=df$Species,featureList = D_list,
                        loss=loss,trace=3,classMethod = 0,stopMethod = 1)

}

##########################################################################
#
#                         REGRESSION TREES
#

rm(list=ls())

data(BostonHousing2)

{
  
# 
# straightLineDistance <- matrix(0,nrow=nrow(BostonHousing2),ncol=nrow(BostonHousing2))
# for(i in 1:(nrow(BostonHousing2)-1))
# {
#   print(i)
#   for(j in (i+1):nrow(BostonHousing2))
#   {
#     d <- distVincentyEllipsoid(BostonHousing2[i,c("lat","lon")],BostonHousing2[j,c("lat","lon")])
#     straightLineDistance[i,j] <- d
#     straightLineDistance[j,i] <- d
#   }
# }


featureList <- with(BostonHousing2,
     list(
       lat=dist(lat),
       lon=dist(lon)
       # spatial = as.dist(straightLineDistance)
     )
)
featureList <- lapply(featureList,as.dist)

value <- BostonHousing2$cmedv

range(BostonHousing2$lat)
range(BostonHousing2$lon)


BostonHousing2[353,"lon"]

library(rpart)

RcppParallel::setThreadOptions(numThreads = 1)

sink(file="log.txt",append = FALSE)


rpart(cmedv~lon+lat,data=BostonHousing2,
      control = rpart.control(cp=0,minbucket = 80,minsplit = 160))

out <- dbrpart(target=value, 
               distanceList  = featureList,
               controls = dbrpart.control(
                              trace=5,
                              minSplit = 160,
                              minBucket = 80,
                              pCenters = 0
                              ))

out2 <- dbrpart(target=value, 
               distanceList  = featureList,
               controls = dbrpart.control(
                 trace=5,
                 minSplit = 160,
                 minBucket = 80,
                 pCenters = 0
               ))


out3 <- dbrpart(target=value, 
                distanceList  = featureList,
                controls = dbrpart.control(
                  trace=5,
                  minSplit = 160,
                  minBucket = 80,
                  pCenters = 0
                ))

sink()

rpart(cmedv~lon+lat,data=BostonHousing2,
      control = rpart.control(cp=0,minbucket = 80,minsplit = 160))->cart

mean(abs(out$y-out$prediction[,1]))
mean(abs(out2$y-out2$prediction[,1]))
mean(abs(out3$y-out3$prediction[,1]))
mean(abs(value-predict(cart)))



out$b
out$nTerm
out$model_risk
out$risk
out$cv_risk

plot(treeGrowth(out))


getTree(out)[,-5] %>% write.table(file="clipboard", sep="\t",row.names = FALSE)




plot(treeGrowth(out),direction="rsquare")



x <- out$classTree
i <- ncol(x)
{
  print(table(x[,i],x[,i-1]))
  i <- i-1
}


# plotcp(cart)

data.frame(variation=c(out$prediction[,1]),
           standard=predict(cart),
           real=value) %>%
  gather("method","predicted",-real) %>%
  ggplot(aes(x=predicted-real,fill=method))+
    geom_density(alpha=0.6)+
    geom_vline(xintercept=0)+
    theme_bw()

mean(abs(value-predict(cart)))

}



##########################################################################
# sink()

# Speed tests here ########################################
RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())
start <- Sys.time()
out <- rcknnpart(featureList = D_list,target=df$Species,loss=loss,trace=1,
                        classMethod = 0,disjointMethod = 0)
end <- Sys.time()
end-start


threshold <- sapply(D_list,function(x){quantile(x[x>0],0.5)})
start <- Sys.time()
out <- rcknnpart(featureList = D_list,target=df$Species,
                        loss=loss,trace=0,
                        disjointMethod = 1,threshold = threshold)
end <- Sys.time()
end-start

RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())
start <- Sys.time()
out <- rcknnpart(featureList = D_list,target=df$Species,loss=loss,trace=0,
                        classMethod = 0,disjointMethod = 3)
end <- Sys.time()
end-start







plot(classTree(out))

varImportance(out)

template <- get_template(out)


png(file="test.png",width=800,height=400)
plot(template,r=0.4,fixRatio = FALSE)
dev.off()

png(file="test2.png",width=800,height=400)
plot(template,r=0.4,fixRatio = TRUE)
dev.off()


# Testing plot_nodes requires some more work
# Would it be a good idea to build a function which
# suggests plotlist?

varImportance(out)

plotList <- list(
  `Petal.Length` = function(id){
    ggplot(iris[as.numeric(id),])+
      geom_boxplot(aes(x=1,y=Petal.Length))+
      theme_bw()+
      scale_y_continuous(limits = c(1,7))
  },
  `Petal.Length:Petal.Width` = function(id)
  {
    ggplot(iris[as.numeric(id),])+
      geom_point(aes(x=Petal.Width,y=Petal.Length))+
      theme_bw()+
      scale_x_continuous(limits = c(1,2.5))+
      scale_y_continuous(limits = c(1,7))
  }
)

target <- function(idIn,idOut)
{
  idIn <- as.numeric(idIn)
  idOut <- as.numeric(idOut)

  df <- iris
  df$id <- 1:nrow(df)
  df <- df[c(idIn,idOut),]

  df$group <- ifelse(df$id %in% idIn,"In","Out")
  ggplot(df,aes(x=group,fill=Species))+
    geom_bar(position="fill",color="black")+
    theme_bw()
}

plots <- plot_nodes(out,target,plotList)

plot(plots[["0"]],rel_widths=c(1,2,1))

plot(template)

pdf(file="test.pdf",width=10,height=4)
plot(
plotTree(template,plots,options=list(tileWidth = 1000,
                                     tileHeight=400,
                                     tileSize = 0.2,
                                     plotTerminal = TRUE))
)
dev.off()

