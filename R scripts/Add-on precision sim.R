
## Simulated precision estimates for TB Vaccine trial add-on study

## Assumptions:
## 13k study population overall including both control and intervention arms
## Annual incidence of cTB will be 0.3%
## Annual incidence of scTB expected to the the same -i.e. ratio is 50/50
## Overall incidence of cTB and scTB is therefore 0.6%

## Incidence of cTB is 0.3%
## Incidence of scTB is 0.3%
## scTB includes progression to cTB, regression to noTB, and maintaining scTB

## Sputum collected 3 or 6 monthly


# Load packages
library(pacman)

p_load(tidyverse,
       lme4,
       broom)


# Set up cohort

cohort <- tibble(id = rep(1:6000, each = 5),
                 month = rep(c(0, 6, 12, 18, 24), times = 6000))


set.seed(180624)

tb_cohort <- cohort |> 
  mutate(scConv = rbinom(n = 30000, size = 1, prob = 0.003))

summarise(tb_cohort, mean(scConv), sum(scConv))

tb_cohort <- tb_cohort |> 
  group_by(id) |> 
  mutate(scTB = sum(scConv))

tb_cohort[tb_cohort$scTB == 1,]








## Simple analysis

# 0.3% chance of cTB at 24 months
# 0.3% chance of scTB at 24 months

simple <- tibble(id = rep(1:6500)) |> 
  mutate(cTB = rbinom(6500, 1, 0.003),
         scTB = rbinom(6500, 1, 0.003),
         arm = case_when(cTB == 1 ~ 1,
                         scTB == 1 ~ 2),
         outcome = case_when(cTB == 1 ~ 1,
                             scTB == 1 ~ 1,
                             .default = 0))

summarise(simple, mean(cTB), mean(scTB), sum(cTB), sum(scTB), sum(outcome))

simple[simple$cTB == 1,]
simple[simple$scTB == 1,]


logit_mod <- glm(outcome ~ as.factor(arm), data = simple, family = binomial(link='logit'))
logit_mod

summary(logit_mod)

exp(logit_mod$coefficients)


log_mod <- glm(outcome ~ as.factor(arm), data = simple, family = binomial(link='identity'))
log_mod

summary(log_mod)










