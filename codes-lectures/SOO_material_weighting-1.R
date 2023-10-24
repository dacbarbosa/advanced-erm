## Chap 2 - Selection on Observable - Weighted regression
## This routine provides examples of manually programmed routines to estimate the
## ATE, and ATT, from weighted regression.
##
## The data are from a version of the Lalonde(1986) data.
##
## At the end of the file, there is an example of the package `weightit`.

## ---- Weighted regression - Manual example with the Lalonde Data

library(dplyr)
library(tidyverse)
library(broom)
library(ggplot2)
library(estimatr)

# For extracting data
library(twang)
data(lalonde)


lalonde$treat_fac <- factor(lalonde$treat,
                            levels = 0:1,
                            labels = c("Untreated", "Treated"))


## ---- Coding it Manually
#An estimate of the propensity score

lreg <- glm(treat ~ age + educ + black + hispan + nodegree +
              married + re74 + re75, data = lalonde, family = binomial(link = 'logit'))
summary(lreg)

p_hat <- predict(lreg, type = 'response')
head(p_hat)

lalonde_augm <- augment(lreg,lalonde, type.predict = 'response')

# Kernel density of propensity score
ggplot(data=lalonde_augm,aes(x = .fitted, colour = treat_fac))+
  geom_density()

# A routine to compute the ATE by IPW (https://www.franciscoyira.com/post/matching-in-r-3-propensity-score-iptw/)
iptw_ate <- function(data,
                     outcome,
                     treatment,
                     prop_score) {
  
  # Renaming the columns for convenience 
  # (i.e. not using {{ }} afterwards)
  data <- data %>% 
    rename(outcome := {{outcome}},
           treatment := {{treatment}},
           prop_score := {{prop_score}})
  
  # Estimation itself
  data %>% 
    mutate(iptw = ifelse(treatment == 1,
                         yes = outcome/prop_score,
                         no = -outcome/(1-prop_score))) %>% 
    pull(iptw) %>% 
    mean()
  
}

# A routine to compute the ATT by IPW (https://www.franciscoyira.com/post/matching-in-r-3-propensity-score-iptw/)
iptw_att <- function(data,
                     outcome,
                     treatment,
                     prop_score) {
  
  # Renaming the columns for convenience 
  # (i.e. not using {{ }} afterwards)
  data <- data %>% 
    rename(outcome := {{outcome}},
           treatment := {{treatment}},
           prop_score := {{prop_score}})
  
  # Estimation itself
  n_treated <- data %>% filter(treatment == 1) %>% nrow()
  
  data %>% 
    mutate(iptw = ifelse(treatment==1,
                         yes = outcome,
                         no = -outcome*prop_score/(1 - prop_score))) %>% 
    pull(iptw) %>% 
    sum() %>% 
    magrittr::divide_by(n_treated)
  
}


# Compute ATE and ATT on the whole sample

iptw_ate(data=lalonde_augm, outcome = re78, treatment = treat, prop_score = .fitted)
iptw_att(data=lalonde_augm, outcome = re78, treatment = treat, prop_score = .fitted)

# Restricting to propensity score between 0.1 and 0.9
lalonde_rest <- lalonde_augm %>%
  filter(.fitted>0.05 & .fitted<0.9)

iptw_ate(data=lalonde_rest, outcome = re78, treatment = treat, prop_score = .fitted)
iptw_att(data=lalonde_rest, outcome = re78, treatment = treat, prop_score = .fitted)



## ----  Normalised IPW (Hirano and Imbens, 2001)
#This estimator normalises the weights based on the sum of the propensity scores
#in the treated and control groups and thus makes the weights themselves
# sum to 1 in each group (so we avoid the awkward situation of weights approaching infinity).

iptw_norm_att <- function(data,
                          outcome,
                          treatment,
                          prop_score) {
  
  # Renaming the columns for convenience 
  # (i.e. not using {{ }} afterwards)
  data <- data %>% 
    rename(outcome := {{outcome}},
           treatment := {{treatment}},
           prop_score := {{prop_score}})
  
  # Estimation itself
  treated_section <- data %>% 
    filter(treatment == 1) %>% 
    summarise(numerator = sum(outcome/prop_score),
              denominator = sum(1/prop_score))
  
  untreated_section <- data %>% 
    filter(treatment == 0) %>% 
    summarise(numerator = sum(outcome/(1-prop_score)),
              denominator = sum(1/(1-prop_score)))
  
  (treated_section$numerator / treated_section$denominator) -
    (untreated_section$numerator / untreated_section$denominator)
}

iptw_norm_att(data=lalonde_augm, outcome = re78, treatment = treat, prop_score = .fitted)




## ---- A routine to bootstrap the result by IPW (https://www.franciscoyira.com/post/matching-in-r-3-propensity-score-iptw/)
library(boot)

# This example is an adaptation of this code script shown in The Effect: https://theeffectbook.net/ch-Matching.html?panelset6=r-code7#panelset6_r-code7 
iptw_boot <- function(data, index = 1:nrow(data)) {
  
  # Slicing with the bootstrap index
  # (this enables the "sampling with replacement")
  data <- data %>%
    slice(index)
  
  # Propensity score estimation
  lreg <- glm(treat ~ age + educ + black + hispan + nodegree +
                married + re74 + re75, data = lalonde, family = binomial(link = 'logit'))
  
  data <- augment(lreg,data, type.predict = 'response')
  
  # Adding the propensity scores to the data
  data <- data %>% 
    rename(prop_score =.fitted) %>% 
    rename(outcome = re78)
  
  # IPTW weighting and estimation
  treated_section <- data %>% 
    filter(treat == 1) %>% 
    summarise(numerator = sum(outcome/prop_score),
              denominator = sum(1/prop_score))
  
  untreated_section <- data %>% 
    filter(treat == 0) %>% 
    summarise(numerator = sum(outcome/(1-prop_score)),
              denominator = sum(1/(1-prop_score)))
  
  (treated_section$numerator / treated_section$denominator) -
    (untreated_section$numerator / untreated_section$denominator)
  
}

b <- boot(lalonde, iptw_boot, R = 100)
b


## - A routine for doubly-robust estimator

# TL:DR, you must change the variables names if you 
# want to re-use this function with your own data

double_robust_estimator <- function(data,
                                    # `index` makes the estimator bootstrap-able
                                    index = 1:nrow(data)) {
  
  data <- data %>% slice(index)
  
  # Getting the "ingredients" to compute the estimator
  # 1. Propensity scores
  # Propensity score estimation
  lreg <- glm(treat ~ age + educ + black + hispan + nodegree +
                married + re74 + re75, data = lalonde, family = binomial(link = 'logit'))
  
  data <- augment(lreg,data, type.predict = 'response')
  
  # Adding the propensity scores to the data
  data <- data %>% 
    rename(prop_score =.fitted) %>% 
    rename(outcome = re78)
  
  # 2. Regression model for treated
  lm_treated <- 
    lm(outcome ~ age + educ + black + hispan + nodegree +
         married + re74 + re75,
       data = data %>% 
         filter(treat == 1))
  
  mu1 <- lm_treated$fitted.values
  
  # 3. Regression model for untreated
  lm_untreated <-
    lm(outcome ~ age + educ + black + hispan + nodegree +
         married + re74 + re75,
       data = data %>%
         filter(treat == 0))
  
  mu0 <- lm_untreated$fitted.values
  
  # 4. Outcomes and Treatment status
  outcome <- data$outcome
  treatment <- data$treat
  prop_score <-- data$prop_score
  
  # 5. Doubly-robust estimator
  estimator <- 
    mean(treatment * (outcome - mu1)/prop_score + mu1) -
    mean((1-treatment)*(outcome - mu0)/(1-prop_score) + mu0)
  
  estimator
}

b_dre <- boot(lalonde, double_robust_estimator, R = 100)

b_dre

## --- Using the package WeightIt
library(cobalt)
library(WeightIt)
#Balancing covariates between treatment groups (binary)
(W.out <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATT"))
summary(W.out)
bal.tab(W.out)

library(survey)
d.w <- svydesign(~1, weights = W.out$weights, data = lalonde)
fit <- svyglm(re78 ~ treat, design = d.w)
summary(fit)
confint(fit)