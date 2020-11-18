# Prototype function for sampling bubble centers
#
# 1. Always start by selecting the record with the largest summed distance from
#    all neighbours. This is the outlieriest option.
#
# 2. Then select the record which is furthest from the one we just selected.
#    
#
# 3. From step 3 onwards, we have two criteria to meet:
#    > Select the record with maximum summed distance from all selected records
#      In essence this selects a ring of outliers
#    > Select the record which is equidistant from all selected records
#      This gives us options in the center as well as the outside 

# One possibility - Convert to separate utility weights (maximise in both cases),
# and then multiply.


# Step 0

set.seed(1)
x <- rnorm(100)
y <- rnorm(100)
col <- rep("black",100)
D <- as.matrix(dist(cbind(x,y)))
plot(x,y,col=col)

# Step 1

sD <- apply(D,1,function(x){sum(sqrt(x))})
col[sD==max(sD)] <- "red"
plot(x,y,col=col)


# Step 2

sD <- D[,col=="red"]
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)


# Step 3
# Complex rules kick in!
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
vD <- 1/apply(D[,col=="red"],1,function(x){var(x)})
pSum <- sD/sum(sD)
pVar <- vD/sum(vD)
p <- pSum+pVar
col[p==max(p[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)


# Step 4
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)


# Step 5
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)



# Step 6
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)



# Step 7
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)



# Step 9
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)


# Step 10
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)

# Step 11
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)


# Step 12
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)



# Step 13
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)




# Step 14
sD <- apply(D[,col=="red"],1,function(x){sum(x)})
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"
plot(x,y,col=col)


##################################################


n<- 200
set.seed(12)

x <- c(rnorm(n),rnorm(n,mean = 0),rnorm(n,mean = 2))
y <- c(rnorm(n),rnorm(n,mean = 2),rnorm(n,mean = 2))

col <- rep("black",length(x))
iter <- rep(-1,length(x))
D <- as.matrix(dist(cbind(x,y)))

# Step 1
sD <- apply(D,1,function(x){sum(sqrt(x))})
col[sD==max(sD)] <- "red"
iter[sD==max(sD)] <- 1

# Step 2
sD <- D[,col=="red"]
iter[sD==max(sD[col=="black"]) & col=="black"] <- 2
col[sD==max(sD[col=="black"]) & col=="black"] <- "red"

for (i in 3:(length(x)-1)){
  sD <- apply(D[,col=="red"],1,function(x){sum(x)})
  vD <- 1/apply(D[,col=="red"],1,function(x){var(x)})
  pSum <- sD/sum(sD)
  pVar <- vD/sum(vD)
  p <- pSum*pVar
  iter[p==max(p[col=="black"]) & col=="black"] <- i
  col[p==max(p[col=="black"]) & col=="black"] <- "red"
}

iter[col=="black"] <- length(x)
col[col=="black"] <- "red"

library(tidyverse)
data.frame(x,y,iter,apply(D,1,sum),col) %>%
  ggplot(aes(x=x,y=y))+
  geom_point(aes(color=iter<100),alpha=0.5)+
  geom_text(aes(label=iter))+
  theme_bw()


data.frame(x,y,iter,d=apply(D,1,sum),col) %>%
  ggplot(aes(x=iter,y=d))+geom_point()
