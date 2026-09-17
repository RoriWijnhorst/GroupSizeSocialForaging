library(rstan)
  
# Load in wrangled data
dfl<- readRDS("Data/Wrangled_data/STAN_input_data.rds")

#Parameters monitored
md <- stan(
  "Models/Mod1.stan",
  data = dfl,
  seed = 123,
  chains = 4,
  iter = 8000,
  warmup = 2000,
  thin = 1,
  cores = 4,
  control = list(adapt_delta = 0.95))

saveRDS(md, "Model_output/Mod1_results.RDS")


