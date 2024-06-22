
### Simulate outcomes from TB Vacc trial with varying levels of specificity of outcome measure

library(pacman)

p_load(tidyverse,
       lme4,
       broom)

## Mimic Franks's simulation

# RR 1.0
# Incidence in BCG arm 0.05
# N=5000 per arm
# NI margin d=1.25 on relative scale

# Number of repetitions and IDs
n_rep <- 1000
n_id <- 10000

set.seed(321)

rr1 <- tibble(
  rep = rep(1:n_rep, each = n_id),
  id = rep(1:n_id, times = n_rep),
  newvac = rep(c(0, 1), each = n_id / 2, times = n_rep),
  rr = 1,
  sens = 1,
  spec = 1,
  spec95 = 0.95,
  case = rbinom(length(case), 1, 0.05)
)

# Calculate obs_case
rr1_obs <- rr1 %>%
  mutate(obs_case = ifelse(
    case == 1,
    rbinom(length(case), 1, rr * sens * 1),    # probability of observing case depends on RR and sens
    rbinom(length(case), 1, (1 - spec))),      # probability of observing no case depends on spec
    obs_case95 = ifelse(
    case == 1,
    rbinom(length(case), 1, rr * sens * 1),    # probability of observing case depends on RR and sens
    rbinom(length(case), 1, (1 - spec95)))     # probability of observing no case depends on spec
    )  

# Validation checks
rr1_obs[4995:5005,]
rr1_obs[9995:10005,]
rr1_obs[14995:15005,]
rr1_obs[19995:110005,]
rr1_obs[114995:115005,]
rr1_obs[809995:8010005,]
rr1_obs[814995:8150005,]

case_rate <- rr1_obs |> 
  group_by(rep, newvac) |> 
  summarise(mean(case), mean(obs_case), mean(obs_case95))

view(case_rate)

rr1_obs |> 
  group_by(newvac) |> 
  summarise(mean(case), mean(obs_case), mean(obs_case95))


## Analysis model
## Log binomial model to compute RR

# Spec = 100%
start.time <- Sys.time()

rr_mod <- rr1_obs |> 
  nest_by(rep) |> 
  mutate(model = list(glm(obs_case ~ as.factor(newvac), data = data, family = binomial(link = "log"))))

rr_mod_ci <- rr_mod |> 
  reframe(tidy(model, conf.int = TRUE, exponentiate = TRUE)) |> 
  filter(term == 'as.factor(newvac)1')

view(rr_mod_ci)

rr_pow <- rr_mod_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high > 1.25 ~ 1,
                           .default = 0))

tab <- table(rr_pow$power)
prop.table(tab)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


# Spec = 95%
rr_mod95 <- rr1_obs |> 
  nest_by(rep) |> 
  mutate(model = list(glm(obs_case95 ~ as.factor(newvac), data = data, family = binomial(link = "log"))))

rr_mod95_ci <- rr_mod95 |> 
  reframe(tidy(model, conf.int = TRUE, exponentiate = TRUE)) |> 
  filter(term == 'as.factor(newvac)1')

view(rr_mod95_ci)

rr_pow95 <- rr_mod95_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high > 1.25 ~ 1,
                           .default = 0))

tab95 <- table(rr_pow95$power)
prop.table(tab95)


# Risk difference: d = 0.024375 as implied by RR of 1.25 ((0.0975*1.25)-0.0975)
# Spec = 95%

start.time <- Sys.time()


rd_mod95 <- rr1_obs |> 
  nest_by(rep) |> 
  mutate(model = list(glm(obs_case95 ~ as.factor(newvac), data = data, family = binomial(link = "identity"))))

rd_mod95_ci <- rd_mod95 |> 
  reframe(tidy(model, conf.int = TRUE, exponentiate = FALSE)) |> 
  filter(term == 'as.factor(newvac)1')

rd_pow95 <- rd_mod95_ci |> 
  select(rep, conf.high, p.value) |> 
  mutate(power = case_when(conf.high > 0.24375 ~ 1,
                           .default = 0))

rd_tab95 <- table(rd_pow95$power)
prop.table(rd_tab95)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken









rr1_obs <- rr1 |> 
  group_by(rep, id) |> 
  mutate(obs_case = case_when(case == 1 ~ rbinom(1, 1, rr*sens*case),       # probability of observing case depends on RR and sens
                              case == 0 ~ rbinom(1, 1, (1-spec)*case))) |>  # probability of observing no case depends on spec
  ungroup()

case_rate <- rr1_obs |> 
  group_by(rep, newvac) |> 
  summarise(mean(case), mean(obs_case))

view(case_rate)







  glm(obs_case ~ as.factor(newvac), data = rr1_obs, family = binomial(link='identity'))
rdiff_mod

rdiff_ci <- tidy(rdiff_mod, conf.int = TRUE, exponentiate = TRUE)
rdiff_ci
