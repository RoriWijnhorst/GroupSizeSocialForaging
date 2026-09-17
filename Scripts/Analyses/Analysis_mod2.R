library(rstan)

# Load in wrangled data
dfl<- readRDS("Data/Wrangled_data/STAN_input_data.rds")

#Parameters monitored
md <- stan(
  "Models/Mod2.stan",
  data = dfl,
  seed = 123,
  chains = 4,
  iter = 8000,
  warmup = 2000,
  thin = 1,
  cores = 4,
  control = list(adapt_delta = 0.95))

md<-"Model_output/Mod2_results.RDS"

