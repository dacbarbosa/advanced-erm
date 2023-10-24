## Chap 2 - Selection on Observable - Estimation of the ATE by nonparametric regression

## This routine simulates non linear potential outcomes model and shows how to compute
## a non parametric regression, and the average treament effect.


library(dplyr)
library(broom)
library(tidyverse)
library(estimatr)
library(ggplot2)

library(stats)
library(locpol)

N<-1000

set.seed(32854516)

## ---- Generate Simulation data
df <- data.frame(x = runif(N))
df$y1 <- df$x^2 - df$x + 1 + rnorm(N, sd = 0.1)
df$y0 <- -0.5*df$x^2 + 1.5*df$x + rnorm(N, sd = 0.5)

df$D <- as.numeric(ifelse(df$x > runif(N) , 1,0))
summary(df$D)

df<-df%>%
  mutate(Y = ifelse(D>0.5,y1,y0))

df<-df%>%
  arrange(x)

# True TE
mean(df$y1 - df$y0)
# Naive TE
mean(df[df$D==1,]$y1) - mean(df[df$D==0,]$y0)

## ---- Fitting a local linear approximation

# Using loess in stats package
r1 <- loess(y1~x,df[df$D==1,],span = 0.75,deg = 1)
m1 <- predict(r1, df$x, se = TRUE,na.omit)
plot(df$x,m1$fit)

# Using locpol package
r1 <- locpol(y1~x,df[df$D==1,],deg = 1,xeval=df$x)
r0 <- locpol(y0~x,df[df$D==0,],deg = 1,xeval=df$x)

summary(r1)

plot(r1)
plot(r0)

m1 <- r1$lpFit$y1
m0 <- r0$lpFit$y0

summary(m0)
summary(m1)

# Append the new moment
df <- df%>%
  arrange(x) %>%
  mutate(m_1 = m1,
         m_0 = m0)

## ---- True Treatment effects
# True TE
mean(df$y1 - df$y0)
# Estimated TE with local polynomial adjustment
mean(df$m_1-df$m_0)

# True TT
mean(df[df$D==1,]$y1-df[df$D==1,]$y0)
# Estimated TT with local polynomial adjustment
mean(df[df$D==1,]$m_1-df[df$D==1,]$m_0)


## ---- Fitting a local linear approximation for the propensity score
prop_hat <- locpol(D~x,data=df,deg = 1,xeval=df$x)
prop_score <- prop_hat$lpFit$D

# Append the new data
df <- df%>%
  arrange(x) %>%
  mutate(ps = prop_score)


## ---- Fitting a local linear approximation for YD and Y(1-D)

df <- df%>%
  mutate(Y1_hat = Y*D,
         Y0_hat = Y*(1-D))

tmpEY1_hat <- locpol(Y1_hat~x,data=df,deg = 1,xeval=df$x)
tmpEY0_hat <- locpol(Y0_hat~x,data=df,deg = 1,xeval=df$x)

# Append the new data
df <- df%>%
  arrange(x) %>%
  mutate(EY1_hat = tmpEY1_hat$lpFit$Y1_hat) %>%
  mutate(EY0_hat = tmpEY0_hat$lpFit$Y0_hat)

df <- df%>%
  filter(ps>0.01 & ps<0.99) %>%
  mutate(mu1_hat = EY1_hat/ps) %>%
  mutate(mu0_hat = EY0_hat/(1-ps))

# True TE
mean(df$y1 - df$y0)
# True TT
mean(df[df$D==1,]$y1-df[df$D==1,]$y0)

# Estimated TE with local polynomial adjustment
mean(df$mu1_hat-df$mu0_hat)
# Estimated TT with local polynomial adjustment
mean(df[df$D==1,]$mu1_hat-df[df$D==1,]$mu0_hat)

# 1- generate a sample with y~X, D: compute the ATT and ATE
# 2- Estimate the ATE and ATT with the local linear regression (EpaK kernel) using a formula similar to m_0 and m_1
# 3- Estimate the ATE and ATT with the local linear regression () using the formula for \mu_0 and \mu_1
# 4- Look at the help for boot(): Write a bootstrap procedure to obtain the s.e. of the estimator. Compare.
# 5- You can do a simulation analysis to compare bs std for both procedures.