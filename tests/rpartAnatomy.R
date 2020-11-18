
library(tidyverse)
library(rpart)
library(rpart.plot)

df <- iris

cartModel <- rpart(Species~.,data=df)

?rpart.object
str(cartModel)
cartModel$frame
cartModel$where
cartModel$call
cartModel$terms
cartModel$splits
cartModel$csplit
cartModel$method
cartModel$cptable
cartModel$variable.importance

summary(cartModel)
