library(cmdstanr)

source("src/scripts/simulatePHECambridge.R")

generatedPHECambridgeData <- generateSimulatedData()
stan_data <- generatedPHECambridgeData$stan_data
stan_data$compute_likelihood = 1
stan_data$debug = 0

#Custom C++ requires cmdstan - Rstan doesn't work in this instance
stanmodel <- cmdstanr::cmdstan_model(
  file.path(getwd(), "src/models/cambridgeModel.stan"),
  include_paths = file.path(getwd(), "src/models"),
  cpp_options = list(
    USER_HEADER = file.path(getwd(), "src/models/dominant_eigenvalue_external.cpp"
  ),
  stanc_options = list("allow-undefined")
)

fit <- stanmodel$sample(
  data = stan_data,
  adapt_delta = 0.99,
  chains = 1, init = 0.1,
  iter_warmup = 100, iter_sampling = 100
)

