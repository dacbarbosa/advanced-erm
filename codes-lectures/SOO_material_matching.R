
## Chap 2: Selection on observable characteristics - Matching

library(dplyr)
library(tidyverse)
library(Matching)
library(cobalt)
library("MatchIt")
library("marginaleffects")

## ---- toy_data_set

set.seed(347628)
n = 1000

## Generating a random data set that complies with a case of selection on observable characteristics

X1 = rnorm(n, 50, 10)
X2 = sample(c(1,0), replace = TRUE, size = n, prob = c(.4,.6)) # binary variable
X3 = rnorm(n, 0, 1) + X1/10
X4 = rnorm(n, 0, 1)

D  = as.numeric(pnorm(rnorm(n, 0, 1) + .01*X1 + .8*X2 + 0.3*X3 + X4 - 2 , 0, 1) > .5)
Y  = rnorm(n, 0, .2) + 1 + 1*D - 0.5*X2 - .4*X3 

# Here the ATE is 

# combining everything in a data frame
df = data.frame(X1, X2, X3, X4, D, Y)

##  ---- Mean by Group and ATE

# Proportion of D=1
count<- table(df$D==1)
perc <- prop.table(count)
cbind(count,perc)

# Compute Mean by group and ATE
df %>%
  group_by(D) %>%
  summarise(mean_Y = mean(Y), number = n() )

# Naive ATE - not the true causal effect
mean(df[df$D==1,]$Y) - mean(df[df$D==0,]$Y)

## ---- Matching the data - using MatchIt Package

# Compute the propensity score

# Using propensity score
match0 <- matchit(D ~ X1 + X2 + X3 + X4, data = df, distance = "glm",
                  method = "nearest", ratio=5, replace = TRUE,
                  exact = ~ X2)
summary(match0)

# Using Mahalanobis distance
match1 <- matchit(D ~ X1 + X2 + X3 + X4, data = df, distance = "mahalanobis",
                  method = "nearest", ratio=5, replace = TRUE,
                  exact = ~ X2)
summary(match1)

## ---- Assessing Match quality

# Balance plot
plot(match1, type = "density", interactive = FALSE,
     which.xs = ~X1 + X3 + X4)
bal.plot(match1, var.name = "X1")

# Graphical representation of mean differnce
love.plot(match0, abs = F)
love.plot(match1, abs = F)

## ---- Estimating the ATT with MatchIt



## Using regression adjustment of Imbens(2015)
# Linear model with covariates  
# calculates the ATT

# Extract the matched data
md <- match.data(match0)
fit1 <- lm(Y ~ D * (X1 + X2 + X3 + X4),
           data = md, weights = weights)

comp1 <- comparisons(fit1,
                     variables = "D",
                     newdata = subset(md, D == 1), # Change here for ATE, or ATNT
                     wts = "weights")

summary(comp1)



## ---- Using blocks on the propensity score estimate

# Subclassification on the PS for the ATT
# requires matching the data on subclasses of the propensity score
mS <- matchit(D ~ X1 + X2 + X3 + X4 , data = df,
              method = "subclass", # Note the new arguments
              estimand = "ATT")

#Extract matched data
md <- match.data(mS)

# Fit with blocks
fitS <- lm(Y ~ subclass * (D * (X1 + X2 + X3 + X4)),
           data = md)

# Compute ATT
compS <- comparisons(fitS,
                     variables = "D",
                     vcov = "HC3",
                     newdata = subset(md, D == 1))
summary(compS)




