# test/run_imifa.R
source(here::here("R", "simfun.R"))

sim <- simfun()

imifa <- run_imifa(
  sim,
  n_iter = 5000,
  burn = 1000,
  thin = 2
)

print(imifa$comparison)