functions {
  //Defined in dominant_eigenvalue_external.cpp
  real dominant_eigenvalue_external(matrix A);

  //Calculate age-group-specific rate of infection
  vector calculateLambda(matrix b, real[] I1, real[] I2){
    int n_age_groups = rows(b);

    vector[n_age_groups] lambda;

    for (i in 1:n_age_groups){
      real lambda_i = 1;
      for (j in 1:n_age_groups){
         lambda_i = lambda_i * (1 - b[i,j])^(I1[j] + I2[j]);
      }
      lambda[i] = 1 - lambda_i;
    }

    return lambda;
  }


  //Integrate SEEII model for single region
  real[,,] simulateOneRegion(real[,] initial_state, matrix[] b, real dI, real dL, real delta_t, int T, int n_age_groups) {

    real state_estimate[6, n_age_groups, T];

    state_estimate[1:5,:,1] = initial_state;

    for (a in 1:n_age_groups){
      state_estimate[6,a,1] = state_estimate[2,a,1] + state_estimate[3,a,1];
    }

    for (t in 2:T){
      vector[n_age_groups] lambda = calculateLambda(b[t], state_estimate[4,:,t-1], state_estimate[5,:,t-1]);

      for (a in 1:n_age_groups){
        state_estimate[6,a,t] = state_estimate[1,a,t-1] * lambda[a] * delta_t; //new infections
        state_estimate[1,a,t] = state_estimate[1,a,t-1] * (1 - lambda[a] * delta_t);
        state_estimate[2,a,t] = state_estimate[2,a,t-1] * (1 - 2 * delta_t/dL) + state_estimate[6,a,t];
        state_estimate[3,a,t] = state_estimate[3,a,t-1] * (1 - 2 * delta_t/dL) + state_estimate[2,a,t-1]*2*delta_t/dL;
        state_estimate[4,a,t] = state_estimate[4,a,t-1] * (1 - 2 * delta_t/dI) + state_estimate[3,a,t-1]*2*delta_t/dL;
        state_estimate[5,a,t] = state_estimate[5,a,t-1] * (1 - 2 * delta_t/dI) + state_estimate[4,a,t-1]*2*delta_t/dI;
      }
    }

    return state_estimate;
  }
}

data {
  int<lower=1> n_beta_pieces;
  int<lower=1> T;
  int<lower=1> n_disease_states;
  int<lower=1> n_regions;
  int<lower=1> n_age_groups;

  int<lower=1, upper = n_age_groups> over75Group; //Index of first group to include over 75s
  int<lower=1> maxL; //Maximum timesteps for infection -> deaths

  int<lower=1> t_lock; //Lockdown time step
  int<lower=1> w_lock; //Lockdown week
  int<lower=1> weeks[T]; //Which week are we in at each time step?

  vector[n_age_groups] p; //Age-group-specific infection-fatality rate
  vector[maxL] pDeathAfterInfection; //Probability of death l days after infection


  real<lower=0> population[n_regions, n_age_groups]; //Per-region and age-group population

  matrix[n_age_groups, n_age_groups] contact_matrix[T]; //Time varying contact matrices

  int<lower=0> deaths[n_regions,n_age_groups,T];

  int<lower=0> n_blood_tests[n_regions, n_age_groups, T];
  int<lower=0> positive_blood_tests[n_regions, n_age_groups, T];

  int<lower=0, upper=1> compute_likelihood;
  int<lower=0, upper=1> debug;
}
transformed data {
  real<lower=0> dL = 4.0;
  real incubationMean = 4.0;
  real incubationSD = 4.0;
  real timeToDeathFromSymptomsMean = 15.0;
  real timeToDeathFromSymptomsSD = 12.1;
  real timeToDeathFromSymptomsShape = timeToDeathFromSymptomsMean*timeToDeathFromSymptomsMean/(timeToDeathFromSymptomsSD*timeToDeathFromSymptomsSD);
  real timeToDeathFromSymptomsRate = timeToDeathFromSymptomsMean/(timeToDeathFromSymptomsSD*timeToDeathFromSymptomsSD);

  real delta_t = 0.5; //Timestep size (days)

  //Need to reverse this to convolve with daily_infections
  vector[maxL] pDeathAfterInfectionReversed = reverse(pDeathAfterInfection);

  matrix[T,T] deathMatrix; //Matrix multiplication  (daily_possible_deaths = deathMatrix * daily_infections) faster than doing in a loop

  {
    row_vector[T] pDeathReversed_row = pDeathAfterInfectionReversed';
    for (t in 1:T){
      int tmpMaxL = maxL;
      if (t <= maxL){
        tmpMaxL = t;

        row_vector[T] tmp;

        tmp[1:t] = pDeathReversed_row[(maxL-t+1):maxL];
        tmp[t+1:T] = rep_row_vector(0.0, T-t);

        deathMatrix[t,] = tmp;
      } else {
        int gap = t - maxL;
        row_vector[T] tmp;
        tmp[1:gap] = rep_row_vector(0.0, gap);
        tmp[(gap+1):(gap+maxL)] = pDeathReversed_row;
        if (t < T){
          tmp[t+1:T] = rep_row_vector(0.0, T-t);
        }
        deathMatrix[t,] = tmp;
      }

    }
  }
  if (debug){
    print(deathMatrix);
  }

}
parameters {
  real<lower=0, upper=1> initial_state_raw[2, n_age_groups, n_regions];

  real<lower=0> dI_raw;

  real<lower=0, upper=1> k_sens; //Serological sensitivity
  real<lower=0, upper=1> k_spec; //Serological specificity

  real<lower=0> sigma_beta; //logscale Step size for random walk of weekly variation sigma_beta
  vector[n_beta_pieces - w_lock] log_beta_steps[n_regions]; //random walk steps

  vector<lower=0>[n_regions] phi_r; //Initial exponential growth per-region
  matrix<lower=0>[n_regions,3] m; //contact matrix multipliers

  real<lower=0> phi_deaths;

}
transformed parameters {
  real initial_state[n_disease_states,n_age_groups,1,n_regions];

  vector[n_beta_pieces] log_beta[n_regions]; //log infection rate
  vector<lower=0>[n_regions] R_0; //Initial R_0 per region

  matrix[n_age_groups, n_age_groups] M[n_regions,2]; //M matrices before and after lockdown (n_regions x 2 array of n_age_groups x n_age_groups matrices)
  matrix[n_age_groups, n_age_groups] C_tilde[n_regions,T]; //C_tilde matrices
  vector<lower=0>[n_regions] Rstar_0; //dominant eigenvalues of initial next-generation matrices
  matrix[n_age_groups, n_age_groups] b[n_regions,T]; //array of infection probability matrices

  real<lower=0> dI;

  vector[T] daily_deaths[n_age_groups, n_regions];
  vector[T] daily_infections[n_age_groups, n_regions];

  real state_estimate[n_disease_states+1, n_age_groups, T, n_regions];

  dI = dI_raw + 2;

  R_0 = phi_r * dI * dot_self(phi_r * dL/2 + 1.0) ./(1 - 1./(dot_self(phi_r * dI/2 + 1)));

  //Initial state taken from Liverpool model for now
  for (a in 1:n_age_groups){
    for (r in 1:n_regions){
      initial_state[1, a, 1, r] = (population[r,a]-5.0)*initial_state_raw[1, a, r] + 1.0;
      initial_state[2, a, 1, r] = (population[r,a]-5.0)*(1.0-initial_state_raw[1, a, r])*initial_state_raw[2, a, r]/2.0 + 1.0;
      initial_state[3, a, 1, r] = (population[r,a]-5.0)*(1.0-initial_state_raw[1, a, r])*initial_state_raw[2, a, r]/2.0 + 1.0;
      initial_state[4, a, 1, r] = (population[r,a]-5.0)*(1.0-initial_state_raw[1, a, r])*(1.0-initial_state_raw[2, a, r])/2.0 + 1.0;
      initial_state[5, a, 1, r] = (population[r,a]-5.0)*(1.0-initial_state_raw[1, a, r])*(1.0-initial_state_raw[2, a, r])/2.0 + 1.0;
    }
  }

  //Random walks on transmission
  for (r in 1:n_regions){
    log_beta[r, 1:w_lock] = rep_vector(0.0, w_lock); //beta_{w_lock,r} = 1
    log_beta[r, (w_lock+1):n_beta_pieces] = cumulative_sum(log_beta_steps[r]);
  }

  //M matrices
  for (r in 1:n_regions){
    //Different susceptibility in over 75s
    //Before Lockdown
    M[r,1][over75Group:n_age_groups,:] = rep_matrix(m[r,1], (n_age_groups-over75Group+1), n_age_groups);
    M[r,1][1:(over75Group-1),:] = rep_matrix(1.0, over75Group-1, n_age_groups);


    //After Lockdown
    M[r,2][over75Group:n_age_groups,:] = rep_matrix(m[r,3], (n_age_groups-over75Group+1), n_age_groups);
    M[r,2][1:(over75Group-1),:] = rep_matrix(m[r,2], over75Group-1, n_age_groups);
  }

  //C_tilde matrices
  for (r in 1:n_regions){
    for (t in 1:(t_lock-1)){
      C_tilde[r,t] = M[r,1] .* contact_matrix[t];
    }
    for (t in t_lock:T){
      C_tilde[r,t] = M[r,2] .* contact_matrix[t];
    }
  }

  //Per-region dominant eigen values of next-generation matrices
  {
    for (r in 1:n_regions){
      matrix[n_age_groups,n_age_groups] Delta; //Next generation matrix for region r
      for (i in 1:n_age_groups){
        Delta[i] = dI * population[r,i] * C_tilde[r,1][i];
      }

      Rstar_0[r] = dominant_eigenvalue_external(Delta);

      if (debug){
        print(Delta);
        print(Rstar_0[r]);
      }
    }
  }

  //Infection probability matrices;
  for (r in 1:n_regions){
    for (t in 1:T){
      b[r,t] = exp(log_beta[r,weeks[t]]) * R_0[r] / Rstar_0[r] * C_tilde[r,t];
    }
  }


  //State estimates from SEEII model
  {
    for (r in 1:n_regions){
      state_estimate[:,:,:,r] = simulateOneRegion(initial_state[:,:,1,r], b[r,:], dI, dL, delta_t, T, n_age_groups);
    }
  }

  //Daily infections stored in state_estimate during integration
  for (r in 1:n_regions){
    for (a in 1:n_age_groups){
      daily_infections[a,r] = to_vector(state_estimate[6,a,:,r]);
    }
  }

  //Daily deaths as lagged proportion of daily infections
  {
    for (r in 1:n_regions){
      for (a in 1:n_age_groups){
        daily_deaths[a,r] = p[a] * deathMatrix * daily_infections[a,r];
      }
    }
  }

  if (debug){
    //print(daily_infections);
    print(daily_deaths);
  }

}

model {

  for (a in 1:n_age_groups){
    for (r in 1:n_regions){
      initial_state_raw[1, a, r] ~ beta(5.0, 0.5);
      initial_state_raw[2, a, r] ~ beta(1.1, 1.1);
    }
  }

  dI_raw ~ gamma(1.43, 0.549);
  phi_r ~ gamma(31.36, 224); //initial exponential growth


  for (j in 1:3){
    for (r in 1:n_regions){
      m[r,j] ~ gamma(4,4); //contact matrix multipliers
    }
  }


  k_sens ~ beta(71.5,29.5); //Serology sensitivity
  k_spec ~ beta(777.5, 9.5); //Serology specificity


  sigma_beta ~ gamma(1,100); //Transmissibility random walk stepsize
  for (r in 1:n_regions){
    log_beta_steps[r,:] ~ normal(0,sigma_beta); //Random walk steps
  }

  //Dispersion parameter for deaths
  phi_deaths ~ gamma(1, 0.2);


  if (compute_likelihood == 1) {
    for (r in 1:n_regions){
      for (a in 1:n_age_groups){

        //Deaths likelihood
        target += neg_binomial_2_lpmf(deaths[r,a] | daily_deaths[a,r], phi_deaths);

        for (t in 1:T){
          //Serum likelihood
          target += binomial_lpmf(positive_blood_tests[r,a,t] | n_blood_tests[r,a,t] , k_sens * (1 - state_estimate[1,a,t,r]/population[r,a]) + (1 - k_spec) * state_estimate[1,a,t,r]/population[r,a]);
        }
      }
    }
  }
}
generated quantities {

}
