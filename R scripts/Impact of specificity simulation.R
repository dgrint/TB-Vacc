
## T1 ==========================================================================
## Simulate outcomes from TB Vacc trial with varying levels of specificity of 
## outcome measure


library(pacman)

p_load(tidyverse,
       lme4,
       broom,
       furrr)


## S1 ==========================================================================
## RR 1.0; True risk 0.02; N=32,000
## NI margin 1.25
## Sensitivity 100%; Specificity 100%


# Number of repetitions and IDs

n_rep <- 2000
n_id <- 32000

# Random number seed

set.seed(102024)

# Simulate data

rr1 <- tibble(
  rep = rep(1:n_rep, each = n_id),
  id = rep(1:n_id, times = n_rep),
  newvac = rep(c(0, 1), each = n_id / 2, times = n_rep),
  rr = 1,
  case = rbinom(n_rep * n_id, 1, 0.02),
  case_n = rbinom(n_rep * n_id, 1, (rr*0.02))
)

rr1 |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(case_n)) 

# Calculate obs_case

rr1_obs <- rr1 %>%
  mutate(obs_case100100 = ifelse(
    case == 1,
    rbinom(length(case), 1, 1),             # probability of observing case depends on RR and sens
    rbinom(length(case), 1, (1 - 1))),      # probability of observing no case depends on spec
    obs_case10095 = ifelse(
      case == 1,
      rbinom(length(case), 1, 1),             # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.95))),   # probability of observing no case depends on spec
    obs_case95100 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 1))),      # probability of observing no case depends on spec
    obs_case9595 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.95))),   # probability of observing no case depends on spec
    obs_case9598 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.98))),    # probability of observing no case depends on spec
    obs_case6485 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.64),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.85)))    # probability of observing no case depends on spec
  ) 

rr1_obs |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(obs_case100100), mean(obs_case10095), mean(obs_case95100), mean(obs_case9595), mean(obs_case9598), 
            mean(obs_case6485)) 

# Validation checks

rr1_obs[4995:5005,]
rr1_obs[9995:10005,]
rr1_obs[14995:15005,]
rr1_obs[19995:110005,]
rr1_obs[114995:115005,]
rr1_obs[809995:8010005,]
rr1_obs[814995:8150005,]


# Use multiple cores

plan(multisession)

## Run analysis models on each rep

## 100% 100%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  # Split reps for parallel processing
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case100100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)



## 100% 95%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case10095 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 100%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case95100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 95%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9595 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)

## 95% 98%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9598 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)



## 64% 85%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  # Split reps for parallel processing
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case6485 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)



## S2 ==========================================================================
## RR 1.25; True risk 0.02; N=32,000
## NI margin 1.25

# Number of repetitions and IDs

n_rep <- 2000
n_id <- 32000

# Random number seed

set.seed(14122024)

# Simulate data

rr125 <- tibble(
  rep = rep(1:n_rep, each = n_id),
  id = rep(1:n_id, times = n_rep),
  newvac = rep(c(0, 1), each = n_id / 2, times = n_rep),
  rr = 1.25,
  case = rbinom(n_rep * n_id, 1, 0.02),
  case_n = rbinom(n_rep * n_id, 1, (rr*0.02))
)

rr125 |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(case_n)) 

# Calculate obs_case

rr125_obs <- rr125 %>%
  mutate(obs_case100100 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 1),           
                    rbinom(length(case_n), 1, (1 - 1))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 1),             
                    rbinom(length(case), 1, (1 - 1))
                  )
           ),
         obs_case10095 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 1),           
                    rbinom(length(case_n), 1, (1 - 0.95))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 1),             
                    rbinom(length(case), 1, (1 - 0.95))
                  )
           ),
         obs_case95100 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 1))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 1))
                  )
           ),
         obs_case9595 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 0.95))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 0.95))
                  )
           ),
         obs_case9598 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 0.98))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 0.98))
                  )
           ),
         obs_case6485 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.64),           
                    rbinom(length(case_n), 1, (1 - 0.85))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.64),             
                    rbinom(length(case), 1, (1 - 0.85))
                  )
           )
  ) 


rr125_obs |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(obs_case100100), mean(obs_case10095), mean(obs_case95100), mean(obs_case9595), mean(obs_case9598), mean(obs_case6485)) 


## Run analysis models on each rep

## 100% 100%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case100100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 100% 95%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case10095 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 100%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case95100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 95%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9595 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 98%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9598 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 64% 85%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case6485 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## S3 ==========================================================================
## RR 1.0; True risk 0.05; N=12,000
## NI margin 1.25
## Sensitivity 100%; Specificity 100%


# Number of repetitions and IDs

n_rep <- 2000
n_id <- 12000

# Random number seed

set.seed(122024)

# Simulate data

rr1 <- tibble(
  rep = rep(1:n_rep, each = n_id),
  id = rep(1:n_id, times = n_rep),
  newvac = rep(c(0, 1), each = n_id / 2, times = n_rep),
  rr = 1,
  case = rbinom(n_rep * n_id, 1, 0.05),
  case_n = rbinom(n_rep * n_id, 1, (rr*0.05))
)

rr1 |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(case_n)) 

# Calculate obs_case

rr1_obs <- rr1 %>%
  mutate(obs_case100100 = ifelse(
    case == 1,
     rbinom(length(case), 1, 1),             # probability of observing case depends on RR and sens
     rbinom(length(case), 1, (1 - 1))),      # probability of observing no case depends on spec
    obs_case10095 = ifelse(
    case == 1,
     rbinom(length(case), 1, 1),             # probability of observing case depends on RR and sens
     rbinom(length(case), 1, (1 - 0.95))),   # probability of observing no case depends on spec
    obs_case95100 = ifelse(
    case == 1,
     rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
     rbinom(length(case), 1, (1 - 1))),      # probability of observing no case depends on spec
    obs_case9595 = ifelse(
    case == 1,
     rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
     rbinom(length(case), 1, (1 - 0.95))),   # probability of observing no case depends on spec
    obs_case9598 = ifelse(
    case == 1,
     rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
     rbinom(length(case), 1, (1 - 0.98))),    # probability of observing no case depends on spec
    obs_case6485 = ifelse(
    case == 1,
     rbinom(length(case), 1, 0.64),          # probability of observing case depends on RR and sens
     rbinom(length(case), 1, (1 - 0.85)))    # probability of observing no case depends on spec
  ) 

rr1_obs |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(obs_case100100), mean(obs_case10095), mean(obs_case95100), mean(obs_case9595), mean(obs_case9598), mean(obs_case6485)) 

# Validation checks

rr1_obs[4995:5005,]
rr1_obs[9995:10005,]
rr1_obs[14995:15005,]
rr1_obs[19995:110005,]
rr1_obs[114995:115005,]
rr1_obs[809995:8010005,]
rr1_obs[814995:8150005,]


# Use multiple cores

plan(multisession)

## Run analysis models on each rep

## 100% 100%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  # Split reps for parallel processing
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case100100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)



## 100% 95%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case10095 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 100%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case95100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 95%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9595 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)

## 95% 98%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9598 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 64% 85%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  # Split reps for parallel processing
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case6485 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)



## S4 ==========================================================================
## RR 1.25; True risk 0.05; N=12,000
## NI margin 1.25

# Number of repetitions and IDs

n_rep <- 2000
n_id <- 12000

# Random number seed

set.seed(1412)

# Simulate data

rr125 <- tibble(
  rep = rep(1:n_rep, each = n_id),
  id = rep(1:n_id, times = n_rep),
  newvac = rep(c(0, 1), each = n_id / 2, times = n_rep),
  rr = 1.25,
  case = rbinom(n_rep * n_id, 1, 0.05),
  case_n = rbinom(n_rep * n_id, 1, (rr*0.05))
)

rr125 |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(case_n)) 

# Calculate obs_case

rr125_obs <- rr125 %>%
  mutate(obs_case100100 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 1),           
                    rbinom(length(case_n), 1, (1 - 1))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 1),             
                    rbinom(length(case), 1, (1 - 1))
                  )
           ),
         obs_case10095 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 1),           
                    rbinom(length(case_n), 1, (1 - 0.95))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 1),             
                    rbinom(length(case), 1, (1 - 0.95))
                  )
           ),
         obs_case95100 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 1))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 1))
                  )
           ),
         obs_case9595 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 0.95))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 0.95))
                  )
           ),
         obs_case9598 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 0.98))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 0.98))
                  )
           ),
         obs_case6485 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.64),           
                    rbinom(length(case_n), 1, (1 - 0.85))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.64),             
                    rbinom(length(case), 1, (1 - 0.85))
                  )
           )
  ) 


rr125_obs |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(obs_case100100), mean(obs_case10095), mean(obs_case95100), mean(obs_case9595), mean(obs_case9598), mean(obs_case6485)) 


## Run analysis models on each rep

## 100% 100%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case100100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 100% 95%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case10095 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 100%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case95100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 95%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9595 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 98%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9598 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 64% 85%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case6485 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## S5 ==========================================================================
## RR 1.0; True risk 0.08; N=7,800
## NI margin 1.25
## Sensitivity 100%; Specificity 100%


# Number of repetitions and IDs

n_rep <- 2000
n_id <- 7300

# Random number seed

set.seed(162024)

# Simulate data

rr1 <- tibble(
  rep = rep(1:n_rep, each = n_id),
  id = rep(1:n_id, times = n_rep),
  newvac = rep(c(0, 1), each = n_id / 2, times = n_rep),
  rr = 1,
  case = rbinom(n_rep * n_id, 1, 0.08),
  case_n = rbinom(n_rep * n_id, 1, (rr*0.08))
)

rr1 |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(case_n)) 

# Calculate obs_case

rr1_obs <- rr1 %>%
  mutate(obs_case100100 = ifelse(
    case == 1,
    rbinom(length(case), 1, 1),             # probability of observing case depends on RR and sens
    rbinom(length(case), 1, (1 - 1))),      # probability of observing no case depends on spec
    obs_case10095 = ifelse(
      case == 1,
      rbinom(length(case), 1, 1),             # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.95))),   # probability of observing no case depends on spec
    obs_case95100 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 1))),      # probability of observing no case depends on spec
    obs_case9595 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.95))),   # probability of observing no case depends on spec
    obs_case9598 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.95),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.98))),    # probability of observing no case depends on spec
    obs_case6485 = ifelse(
      case == 1,
      rbinom(length(case), 1, 0.64),          # probability of observing case depends on RR and sens
      rbinom(length(case), 1, (1 - 0.85)))    # probability of observing no case depends on spec
  ) 

rr1_obs |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(obs_case100100), mean(obs_case10095), mean(obs_case95100), mean(obs_case9595), mean(obs_case9598), mean(obs_case6485)) 

# Validation checks

rr1_obs[4995:5005,]
rr1_obs[9995:10005,]
rr1_obs[14995:15005,]
rr1_obs[19995:110005,]
rr1_obs[114995:115005,]
rr1_obs[809995:8010005,]
rr1_obs[814995:8150005,]


# Use multiple cores

plan(multisession)

## Run analysis models on each rep

## 100% 100%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  # Split reps for parallel processing
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case100100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)



## 100% 95%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case10095 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 100%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case95100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 95%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9595 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)

## 95% 98%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9598 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 64% 85%

start.time <- Sys.time()

rr_mod_ci <- rr1_obs |> 
  # Split reps for parallel processing
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case6485 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Power

tab <- table(rr_pow$power)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)



## S6 ==========================================================================
## RR 1.25; True risk 0.08; N=7,300
## NI margin 1.25

# Number of repetitions and IDs

n_rep <- 2000
n_id <- 7300

# Random number seed

set.seed(1812)

# Simulate data

rr125 <- tibble(
  rep = rep(1:n_rep, each = n_id),
  id = rep(1:n_id, times = n_rep),
  newvac = rep(c(0, 1), each = n_id / 2, times = n_rep),
  rr = 1.25,
  case = rbinom(n_rep * n_id, 1, 0.08),
  case_n = rbinom(n_rep * n_id, 1, (rr*0.08))
)

rr125 |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(case_n)) 

# Calculate obs_case

rr125_obs <- rr125 %>%
  mutate(obs_case100100 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 1),           
                    rbinom(length(case_n), 1, (1 - 1))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 1),             
                    rbinom(length(case), 1, (1 - 1))
                  )
           ),
         obs_case10095 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 1),           
                    rbinom(length(case_n), 1, (1 - 0.95))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 1),             
                    rbinom(length(case), 1, (1 - 0.95))
                  )
           ),
         obs_case95100 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 1))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 1))
                  )
           ),
         obs_case9595 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 0.95))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 0.95))
                  )
           ),
         obs_case9598 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.95),           
                    rbinom(length(case_n), 1, (1 - 0.98))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.95),             
                    rbinom(length(case), 1, (1 - 0.98))
                  )
           ),
         obs_case6485 = 
           ifelse(newvac == 1,
                  ifelse(
                    case_n == 1,
                    rbinom(length(case_n), 1, 0.64),           
                    rbinom(length(case_n), 1, (1 - 0.85))
                  ),
                  ifelse(
                    case == 1,
                    rbinom(length(case), 1, 0.64),             
                    rbinom(length(case), 1, (1 - 0.85))
                  )
           )
  ) 


rr125_obs |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(obs_case100100), mean(obs_case10095), mean(obs_case95100), mean(obs_case9595), mean(obs_case9598), mean(obs_case6485)) 


## Run analysis models on each rep

## 100% 100%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case100100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 100% 95%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case10095 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 100%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case95100 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 95%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9595 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 95% 98%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case9598 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)


## 64% 85%

start.time <- Sys.time()

rr_mod_ci <- rr125_obs |> 
  group_split(rep) |> 
  future_map(~ {
    model <- glm(obs_case6485 ~ as.factor(newvac), data = ., family = binomial(link = "log"))
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  }) |> 
  # Combine results into one data frame
  bind_rows(.id = "rep") |> 
  filter(term == 'as.factor(newvac)1')

end.time <- Sys.time()
end.time - start.time


rr_alpha <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(alpha = case_when(conf.high <= 1.25 ~ 1,
                           .default = 0))

# Alpha

tab <- table(rr_alpha$alpha)
prop <- prop.table(tab)
prop

# Monte Carlo error

sqrt((prop[2]*(prop[1]))/n_rep)
















