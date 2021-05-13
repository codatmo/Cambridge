#Calculate age-group-specific rate of infection
calculateLambda <- function(b, I1, I2){
  nAgeGroups <- nrow(b)
  lambda <- array(0,dim=c(nAgeGroups))
  for (i in 1:nAgeGroups){
    lambda_i <- 1
    for (j in 1:nAgeGroups){
      lambda_i <- lambda_i * (1 - b[i,j]) ^ (I1[j] + I2[j])
    }
    lambda[i] = 1 - lambda_i
  }
  lambda
}


#Function integrating SEEII model for a single region
simulateOneRegion <- function(
  initial_state, b, dI, dL, delta_t, maxTime, nAgeGroups
){
  state_estimate <- array(NaN, dim = c(6, nAgeGroups, maxTime))

  state_estimate[1:5, , 1] = initial_state


  #Daily infections at initial state is initial number of infections?
  for (a in 1:nAgeGroups){
    state_estimate[6,a,1] = state_estimate[2,a,1] + state_estimate[3,a,1]
  }

  for (t in 2:maxTime){
    lambda = calculateLambda(b[,,t], state_estimate[4,,t-1], state_estimate[5,,t-1])

    for (a in 1:nAgeGroups){
      state_estimate[6,a,t] = state_estimate[1,a,t-1] * lambda[a] * delta_t
      state_estimate[1,a,t] = state_estimate[1,a,t-1] * (1 - lambda[a] * delta_t)
      state_estimate[2,a,t] = state_estimate[2,a,t-1] * (1 - 2 * delta_t/dL) + state_estimate[6,a,t]
      state_estimate[3,a,t] = state_estimate[3,a,t-1] * (1 - 2 * delta_t/dL) + state_estimate[2,a,t-1]*2*delta_t/dL
      state_estimate[4,a,t] = state_estimate[4,a,t-1] * (1 - 2 * delta_t/dI) + state_estimate[3,a,t-1]*2*delta_t/dL
      state_estimate[5,a,t] = state_estimate[5,a,t-1] * (1 - 2 * delta_t/dI) + state_estimate[4,a,t-1]*2*delta_t/dI
    }
  }

  state_estimate
}


#Simulate from prior for PHE/Cambridge model
generateSimulatedData <- function(
  seed = 0,
  nRegions = 2,
  nAgeGroups = 3,
  over75Group = 3, #Which age group is the first of the over 75s
  population = 1000000 #Population of each age group in each region (not realistic)
){

  set.seed(seed)
  initial_time <- 0
  n_disease_states <- 5

  delta_t <- 0.5 #Time steps are half-days
  nBetaPieces <- 20 #20 weeks
  maxTime = nBetaPieces*14 #Number of half-days
  wLock <- 6 #Lockdown starting at 6th week
  tLock <- wLock*14+1 #Time step of lockdown start
  weeks <- rep(seq(1,nBetaPieces), each = 14) #Indicator for which week each half-day belongs to

  #In reality these should come from e.g. google mobility data
  contactMatrices <- array(1, dim = c(nAgeGroups, nAgeGroups, maxTime)) * 1e-7

  population <- array(population, dim = c(nRegions, nAgeGroups))

  #Assume 1000 blood tests in each region and age group
  n_blood_tests <- array(1000, dim = c(nRegions, nAgeGroups, maxTime))


  ifr <- c(rep(c(0.001), each = nAgeGroups-1), 0.05) #Low ifr except for over 75s for example

  #Prior parameters for dL, dI, dT
  mu_dL = 4.00;
  sigma_dL = 0.20;
  mu_dI = 3.06;
  sigma_dI = 0.21;
  mu_dT = 16.00;
  sigma_dT = 0.71;

  timeToDeathFromSymptomsMean = 15.0
  timeToDeathFromSymptomsSD = 12.1
  timeToDeathFromSymptomsShape = timeToDeathFromSymptomsMean*timeToDeathFromSymptomsMean/(timeToDeathFromSymptomsSD*timeToDeathFromSymptomsSD)
  timeToDeathFromSymptomsRate = timeToDeathFromSymptomsMean/(timeToDeathFromSymptomsSD*timeToDeathFromSymptomsSD)

  maxL <- 100

  pDeathAfterInfection <- sapply(seq(0,maxL-2,by=0.5), function(d){
    pgamma(d+0.5,timeToDeathFromSymptomsShape, timeToDeathFromSymptomsRate)-pgamma(d,timeToDeathFromSymptomsShape, timeToDeathFromSymptomsRate)
  })
  pDeathAfterInfection <- c(pDeathAfterInfection, 1 - sum(pDeathAfterInfection))
  maxL <- length(pDeathAfterInfection)

  #Draw from priors from model
  initial_state_raw <- array(0, dim = c(2, nAgeGroups, nRegions))
  for (a in seq(1,nAgeGroups)){
    for (r in seq(1,nRegions)){
      initial_state_raw[,a,r] <- c(rbeta(1, 5.0, 0.5), rbeta(1, 1.1, 1.1))
    }
  }


  dI <- rgamma(1,1.43, 0.549) + 2
  phi_r <- rgamma(nRegions, 31.36, 224) #initial exponential growth

  #Contact Matrix Multipliers
  contactMatrixMultipliers <- array(0, dim = c(nRegions,3))
  for (r in seq(1,nRegions)){
    contactMatrixMultipliers[r,] <- rgamma(3,4,4)
  }

  #Transmissibility random walk step size
  sigma_beta <- rgamma(1,1,100)
  log_beta_steps <- matrix(rnorm(nRegions*(nBetaPieces-wLock), 0, sigma_beta), nrow = nRegions)

  #Deaths dispersion
  phi_deaths <- rgamma(1, 1, 0.2)

  #Blood test sensitivity/specificity
  k_sens <- rbeta(1, 71.5, 29.5) #Serology sensitivity
  k_spec <- rbeta(1, 777.5, 9.5) #Serology specificity

  dL <-  rnorm(1, mu_dL, sigma_dL)
  dT <- rnorm(1, mu_dT, sigma_dT)

  #R_0
  R_0 <- phi_r * dI * (phi_r * dL/2 + 1)^2 / (1 - 1/(phi_r * dI/2 + 1)^2)

  #Initial state
  initial_state <- array(0, dim = c(5,nAgeGroups,1,nRegions))
  for (a in 1:nAgeGroups){
    for (r in 1:nRegions){
      initial_state[1,a,1,r] <- (population[r,a] - 5.0) * initial_state_raw[1,a,r] + 1.0
      initial_state[2,a,1,r] <- (population[r,a] - 5.0) * (1 - initial_state_raw[1,a,r]) * initial_state_raw[2,a,r]/2.0 + 1.0
      initial_state[3,a,1,r] <- (population[r,a] - 5.0) * (1 - initial_state_raw[1,a,r]) * initial_state_raw[2,a,r]/2.0 + 1.0
      initial_state[4,a,1,r] <- (population[r,a] - 5.0) * (1 - initial_state_raw[1,a,r]) * (1 - initial_state_raw[2,a,r])/2.0 + 1.0
      initial_state[5,a,1,r] <- (population[r,a] - 5.0) * (1 - initial_state_raw[1,a,r]) * (1 - initial_state_raw[2,a,r])/2.0 + 1.0
    }
  }


  #Log infection rate
  log_beta <- array(NaN, dim = c(nRegions,nBetaPieces))
  for (r in 1:nRegions){
    log_beta[r, 1:wLock] <- 0.0
    log_beta[r, (wLock++1):nBetaPieces] <- cumsum(log_beta_steps[r,])
  }


  #M matrices
  M <- array(NaN, dim = c(nAgeGroups, nAgeGroups, nRegions, 2))
  for (r in 1:nRegions){
    #Different susceptibilities in over 75s
    #Before lockdown
    M[over75Group:nAgeGroups,,r,1] <- array(contactMatrixMultipliers[r,1], dim = c(nAgeGroups-over75Group+1, nAgeGroups))
    M[1:(over75Group-1),,r,1] <- array(contactMatrixMultipliers[r,2], dim = c(over75Group-1, nAgeGroups))

    #Before lockdown
    M[over75Group:nAgeGroups,,r,2] <- array(contactMatrixMultipliers[r,3], dim = c(nAgeGroups-over75Group+1, nAgeGroups))
    M[1:(over75Group-1),,r,2] <- array(contactMatrixMultipliers[r,2], dim = c(over75Group-1, nAgeGroups))

  }


  C_tilde <- array(NaN, dim = c(nAgeGroups, nAgeGroups, nRegions, maxTime))
  for (r in 1:nRegions){
    for (t in 1:(tLock-1)){
      C_tilde[,,r,t] = M[,,r,1] * contactMatrices[,,t]
    }

    for (t in tLock:maxTime){
      C_tilde[,,r,t] = M[,,r,2] * contactMatrices[,,t]
    }
  }

  #Initial R
  Rstar_0 <- array(0, dim = c(nRegions))
  for (r in 1:nRegions){
    Delta <- array(0, dim = c(nAgeGroups, nAgeGroups))
    for(i in 1:nAgeGroups){
      Delta[i,] <- dI * population[r,i] * C_tilde[i,,r,1]
    }

    Rstar_0[r] <- Re(eigen(Delta)$values[1])
  }

  #Infection probability matrices
  b <- array(0, dim = c(nAgeGroups, nAgeGroups, nRegions, maxTime))
  for (r in 1:nRegions){
    for (t in 1:maxTime){
      b[,,r,t] <- exp(log_beta[r,weeks[t]]) * R_0[r] / Rstar_0[r] * C_tilde[,,r,t]
    }
  }

  state_estimate <- array(NaN, dim = c(6,nAgeGroups,maxTime,nRegions))
  for (r in 1:nRegions){
    state_estimate[,,,r] <- simulateOneRegion(initial_state[,,,r], b[,,r,], dI, dL, delta_t, maxTime, nAgeGroups)
  }

  daily_infections <- state_estimate[6,,,]

  #Daily deaths as lagged proportion of infections
  daily_deaths <- array(0,dim = c(nAgeGroups, maxTime, nRegions))
  for (r in 1:nRegions){
    for (a in 1:nAgeGroups){
      for (t in 1:maxTime){

        cumulativePossibleDeaths = 0
        tmpMaxL = maxL
        if (t < tmpMaxL){
          tmpMaxL = t
        }
        cumulativePossibleDeaths <- sum(rev(pDeathAfterInfection)[(maxL-tmpMaxL+1):maxL] * daily_infections[a,seq(t-tmpMaxL+1,t),r])

        daily_deaths[a,t,r] = ifr[a] * cumulativePossibleDeaths
      }
    }
  }

  deaths <- apply(daily_deaths, c(1,2,3), function(x) { rnbinom(1, size = phi_deaths, mu = x)})

  positive_blood_tests <- array(0,dim = c(nRegions, nAgeGroups, maxTime))
  for (r in 1:nRegions){
    for (a in 1:nAgeGroups){
      for (t in 1:maxTime){
        positive_blood_tests[r,a,t] <- rbinom(1, n_blood_tests[r,a,t] , k_sens * (1 - state_estimate[1,a,t,r]/population[r,a]) + (1 - k_spec) * state_estimate[1,a,t,r]/population[r,a])
      }
    }
  }

  contactMatrices_reshaped <- array(0, dim = c(dim(contactMatrices)[3],nAgeGroups, nAgeGroups))
  for (i in 1:dim(contactMatrices)[3]){
    contactMatrices_reshaped[i,,] <- contactMatrices[,,i]
  }

  #Data for Stan model
  stan_data = list(

    n_beta_pieces = nBetaPieces,
    T = maxTime,
    n_disease_states = 5,
    n_regions = nRegions,
    n_age_groups = nAgeGroups,
    over75Group = over75Group,

    maxL = length(pDeathAfterInfection),

    t_lock = tLock,
    w_lock = wLock,
    weeks = weeks,

    p = ifr, #Age-group specific ifr
    pDeathAfterInfection = pDeathAfterInfection, #f[i] = Probability of death i days after infection

    population = population, #Agre and region populations
    contact_matrix = contactMatrices_reshaped,

    deaths = aperm(deaths,c(3,1,2)),

    n_blood_tests = n_blood_tests,
    positive_blood_tests = positive_blood_tests
  )

  #Return the ground truth values to compare against
  ground_truth = list(
    state_estimate = state_estimate,
    daily_infections = daily_infections,
    daily_deaths = daily_deaths,
    log_beta = log_beta,
    b = b,
    C_tilde = C_tilde,
    contactMatrixMultipliers = contactMatrixMultipliers,
    ifr = ifr,
    pDeathAfterInfection = pDeathAfterInfection
  )

  list(stan_data = stan_data, ground_truth = ground_truth)

}
