data {
  int<lower=1> J;
  int<lower=1> N_z1;
  array[N_z1] int<lower=1> ind_z1;
  int<lower=1> N_z2;
  array[N_z2] int<lower=1> ind_z2;

  int<lower=1> Ng;
  int<lower=1> Nt_t;
  int<lower=1> Nt_g;

  vector[N_z1] Day;
  vector[N_z1] Seq_t;
  vector[N_z1] Grp_t;
  vector[N_z2] Seq_g;
  vector[N_z2] Grp_g;

  array[N_z1] int<lower=0> Scr_opp1;
  array[N_z1] int<lower=0> For_opp1;
  array[N_z2] int<lower=0> Scr_opp2;
  array[N_z2] int<lower=0> For_opp2;

  array[N_z1] int<lower=1> Group_t;
  array[N_z1] int<lower=1> Trial_t;
  array[N_z2] int<lower=1> Group_g;
  array[N_z2] int<lower=1> Trial_g;

  array[N_z1] int<lower=0> Scr1;
  array[N_z1] int<lower=0> n_events1;
  array[N_z2] int<lower=0> Scr2;
  array[N_z2] int<lower=0> n_events2;
}

transformed data {
  vector[N_z1] Prop_opp_obs1;
  vector[N_z2] Prop_opp_obs2;

  vector[N_z1] Prop_opp_c1;
  vector[N_z2] Prop_opp_c2;

  for(i in 1:N_z1){
    Prop_opp_obs1[i] = Scr_opp1[i] * 1.0 / For_opp1[i];
  }

  for(i in 1:N_z2){
    Prop_opp_obs2[i] = Scr_opp2[i] * 1.0 / For_opp2[i];
  }

  Prop_opp_c1 = Prop_opp_obs1 - mean(Prop_opp_obs1);
  Prop_opp_c2 = Prop_opp_obs2 - mean(Prop_opp_obs2);
}

parameters {
  real mu0_1;
  real mu0_2;

  vector[7] beta_t;
  vector[3] beta_g;

  real beta_p1;
  real beta_p2;

  vector<lower=0>[2] sd_int;
  matrix[2, J] z_int;
  cholesky_factor_corr[2] R_int_chol;

  vector[Nt_t] zTrial_t;
  real<lower=0> sigma_Trial_t;
  
  vector[Ng] zGroup_t;
  real<lower=0> sigma_Group_t;

  vector[Nt_g] zTrial_g;
  real<lower=0> sigma_Trial_g;
  vector[Ng] zGroup_g;
  real<lower=0> sigma_Group_g;

}

transformed parameters {
  matrix[J,2] Int_RE =
  (diag_pre_multiply(sd_int, R_int_chol) * z_int)';

  vector[Nt_t] tr_t = zTrial_t * sigma_Trial_t;
  vector[Ng]   gr_t = zGroup_t * sigma_Group_t;

  vector[Nt_g] tr_g = zTrial_g * sigma_Trial_g;
  vector[Ng]   gr_g = zGroup_g * sigma_Group_g;
}

model {
  vector[J] a1 = col(Int_RE, 1);
  vector[J] a2 = col(Int_RE, 2);

  vector[N_z1] mu1;
  vector[N_z2] mu2;

  mu1 = mu0_1 + a1[ind_z1]
       + beta_t[1] * Day
       + beta_t[2] * Seq_t
       + beta_t[3] * Grp_t
       + beta_t[4] * Day .* Seq_t
       + beta_t[5] * Day .* Grp_t
       + beta_t[6] * Seq_t .* Grp_t
       + beta_t[7] * Day .* Seq_t .* Grp_t
       + beta_p1.* Prop_opp_c1
       + tr_t[Trial_t]
       + gr_t[Group_t];

  mu2 = mu0_2 + a2[ind_z2]
       + beta_g[1] * Seq_g
       + beta_g[2] * Grp_g
       + beta_g[3] * Seq_g .* Grp_g
       + beta_p2 .* Prop_opp_c2
       + tr_g[Trial_g]
       + gr_g[Group_g];

  Scr1 ~ binomial_logit(n_events1, mu1);
  Scr2 ~ binomial_logit(n_events2, mu2);

  mu0_1 ~ normal(0, 1);
  mu0_2 ~ normal(0, 1);
  beta_t ~ normal(0, 1);
  beta_g ~ normal(0, 1);
  beta_p1 ~ normal(0, 1);
  beta_p2 ~ normal(0, 1);

  to_vector(zTrial_t) ~ std_normal();
  sigma_Trial_t ~ exponential(2);

  to_vector(zGroup_t) ~ std_normal();
  sigma_Group_t ~ exponential(2);

  to_vector(zTrial_g) ~ std_normal();
  sigma_Trial_g ~ exponential(2);

  to_vector(zGroup_g) ~ std_normal();
  sigma_Group_g ~ exponential(2);

  sd_int ~ exponential(2);
  R_int_chol ~ lkj_corr_cholesky(2);
  to_vector(z_int) ~ std_normal();
}

generated quantities {
  matrix[2,2] R_int = R_int_chol * R_int_chol';
  matrix[2,2] S_int = diag_matrix(sd_int);
  matrix[2,2] P_int = S_int * R_int * S_int;
  
  vector[N_z1] x_obs1;
  vector[N_z2] x_obs2;
  
  for(i in 1:N_z1){
    x_obs1[i] = Scr_opp1[i] * 1.0 / For_opp1[i];
  }

  for(i in 1:N_z2){
    x_obs2[i] = Scr_opp2[i] * 1.0 / For_opp2[i];
  }
  
  real V_int1   = P_int[1,1];
  real V_int2   = P_int[2,2];
  real V_trial1 = sigma_Trial_t^2; 
  real V_trial2 = sigma_Trial_g^2; 
  real V_group1 = sigma_Group_t^2; 
  real V_group2 = sigma_Group_g^2; 

  real cor_int1_int2     = R_int[1,2];

  real V_dist = pi() * pi() / 3;

  
vector[N_z1] eta_social1;
  for (i in 1:N_z1)
  eta_social1[i] = beta_p1 * Prop_opp_c1[i];

real V_social1 = variance(eta_social1);

vector[N_z2] eta_social2; 
  for (i in 1:N_z2) { eta_social2[i] = beta_p2 * Prop_opp_c2[i]; } 
real V_social2 = variance(eta_social2); 


// Total variance (design fixed effects removed)
real V_total1 = V_int1 + V_social1 + V_trial1 + V_group1 + V_dist;
real V_total2 = V_int2 + V_social2 + V_trial2 + V_group2 + V_dist;


// ---------------------------
// Repeatabilities
// ---------------------------

real VC_int1   = V_int1   / V_total1;
real VC_int2   = V_int2   / V_total2;

real VC_trial1 = V_trial1 / V_total1;
real VC_trial2 = V_trial2 / V_total2;

real VC_group1 = V_group1 / V_total1;
real VC_group2 = V_group2 / V_total2;

real VC_social1 = V_social1 / V_total1;
real VC_social2 = V_social2 / V_total2;

real VC_residual1 = V_dist /  V_total1;
real VC_residual2 = V_dist /  V_total2;

real VC_within1 = (V_social1+V_dist) /  V_total1;
real VC_within2 = (V_social2+V_dist) /  V_total2;
}
