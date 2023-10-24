
## Chap 2: Selection on observable characteristics - Linear Model

## This routine simulates simple data with selection on observables.
## It shows how to conduct estimation of the treatment effect parameters using
## a linear model.

library(dplyr)
library(tidyverse)
library(estimatr)
library(ggplot2)


## ---- three_datasets
set.seed(347628)
n = 10000

#Simulate three datasets:
# Common observables
X1 = rnorm(n, 50, 10)
X2 = sample(c(1,0), replace = TRUE, size = n, prob = c(.4,.6)) # binary variable
X3 = rnorm(n, 0, 1) + X1/10
u_Y =  rnorm(n, 0, .2)
u_D = rnorm(n, 0, 1)
  
# one with selection on observable, but observable uncorrelated with the outcome
# Hence it does not matter for the ATE

D  = as.numeric(pnorm(u_D + .01*X1 - 1 , 0, 1) > .5)
# Same as: D  = as.numeric(rnorm(n, 0, 1) + .01*X1 - 1 , 0, 1) > 0)
Y  = u_Y + 1 + 1*D - 0.5*X2  

# combining everything in a data frame
df_ch2_data1 = data.frame(X1, X2, X3, D, Y)


# one with selection on observable, and correlated observable
D  = as.numeric(pnorm(u_D + .01*X1 + .8*X2 - 1 , 0, 1) > .5)
Y  = u_Y + 1 + 1*D - 0.5*X2  

# combining everything in a data frame
df_ch2_data2 = data.frame(X1, X2, X3, D, Y)

# one with worse selection on observable, correlated observable.
D  = as.numeric(pnorm(u_D + .01*X1 + .8*X2 - 1 , 0, 1) > .5)
Y  = u_Y + 1 + 1*D - 0.5*X2 - .4*X1 

# combining everything in a data frame
df_ch2_data3 = data.frame(X1, X2, X3, D, Y)

# one with selection on observable, and no overlap

D  = as.numeric(pnorm(u_D + .01*X1 - 1 , 0, 1) > .5)
D = 0*as.numeric(X1>50) + D *as.numeric(X1<=50)
Y  = u_Y + 1 + 1*D - 0.5*X2  


# combining everything in a data frame
df_ch2_data4 = data.frame(X1, X2, X3, D, Y)

ggplot(data = df_ch2_data4) +
 geom_point(mapping= aes(x = X1, y = D))
  

# Here the ATE is 1

##  ---- group_mean

#Mean by Group and ATE

# Proportion of D=1
count<- table(df_ch2_data1$D==1)
perc <- prop.table(count)
cbind(count,perc)

# Compute Mean by group and ATE
df_ch2_data1 %>%
  group_by(D) %>%
  summarise(mean_Y = mean(Y), number = n() )

# Naive ATE - not the true causal effect
mean(df_ch2_data1[df_ch2_data1$D==1,]$Y) - mean(df_ch2_data1[df_ch2_data1$D==0,]$Y)
mean(df_ch2_data2[df_ch2_data2$D==1,]$Y) - mean(df_ch2_data2[df_ch2_data2$D==0,]$Y)
mean(df_ch2_data1[df_ch2_data3$D==1,]$Y) - mean(df_ch2_data1[df_ch2_data3$D==0,]$Y)
mean(df_ch2_data2[df_ch2_data4$D==1,]$Y) - mean(df_ch2_data2[df_ch2_data4$D==0,]$Y)



## ---- Linear_Model

# - Data Frame 2
# Centering the variables
demeaned_df_ch2_data2<- df_ch2_data2 %>%
  mutate(x1 = X1 - mean(X1),x2 = X2 - mean(X2),x3 = X3 - mean(X3))

reg_pooled_model <- lm_robust(Y ~ D+ x1 + x2 + x3, demeaned_df_ch2_data2)
summary(reg_pooled_model)

# - Data frame 3
# Centering the variables
demeaned_df_ch2_data3<- df_ch2_data3 %>%
  mutate(x1 = X1 - mean(X1),x2 = X2 - mean(X2),x3 = X3 - mean(X3))

reg_pooled_model <- lm(Y ~ D+ x1 + x2 + x3, demeaned_df_ch2_data3)
summary(reg_pooled_model)$ coefficients


