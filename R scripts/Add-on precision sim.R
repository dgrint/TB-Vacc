
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

cTB <- tibble(id = rep(1:6500)) |> 
  mutate(arm = 'cTB',
         cTB = rbinom(6500, 1, 0.003),
         outcome = case_when(cTB == 1 ~ 1,
                             cTB == 0 ~ 0))

summarise(cTB, mean(cTB), sum(cTB), sum(outcome))

scTB <- tibble(id = rep(1:6500)) |> 
  mutate(arm = 'scTB',
         scTB = rbinom(6500, 1, 0.003),
         outcome = case_when(scTB == 1 ~ 1,
                             scTB == 0 ~ 0))

summarise(scTB, mean(scTB), sum(scTB), sum(outcome))

# Bind together
simple <- bind_rows(cTB, scTB)



# Risk difference model
rdiff_mod <- glm(outcome ~ as.factor(arm), data = simple, family = binomial(link='identity'))
rdiff_mod

rdiff_ci <- tidy(rdiff_mod, conf.int = TRUE)
rdiff_ci

# Prevalence ratio model
log_mod <- glm(outcome ~ as.factor(arm), data = simple, family = binomial(link='log'))
log_mod

log_ci <- tidy(log_mod, conf.int = TRUE, exponentiate = TRUE)
log_ci





# 0.3% chance of scTB at 12 months
# 50% chance of scTB -> cTB at 24 months

m12 <- tibble(id = rep(1:6500)) |> 
  mutate(arm = 'scTB',
         scTB = rbinom(6500, 1, 0.003))

m24 <- m12 |> 
  filter(scTB == 1 ) |> 
  group_by(id) |> 
  mutate(cTB = rbinom(1, 1, 0.5),
         outcome = case_when(cTB == 1 ~ 1,
                             cTB == 0 ~ 0)) |> 
  ungroup()

get_wilson_CI <- function(x, alpha = 0.05) {
  #-----------------------------------------------------------------------------
  # Compute the Wilson (aka Score) confidence interval for a popn. proportion
  #-----------------------------------------------------------------------------
  # x        vector of data (zeros and ones)
  # alpha    1 - (confidence level)
  #-----------------------------------------------------------------------------
  n <- length(x)
  p_hat <- mean(x)
  SE_hat_sq <- p_hat * (1 - p_hat) / n
  crit <- qnorm(1 - alpha / 2)
  omega <- n / (n + crit^2)
  A <- p_hat + crit^2 / (2 * n)
  B <- crit * sqrt(SE_hat_sq + crit^2 / (4 * n^2))
  CI <- c('lower' = omega * (A - B), 
          'upper' = omega * (A + B))
  return(CI)
}

summarise(m24, mean(outcome))
get_wilson_CI(m24$outcome)


