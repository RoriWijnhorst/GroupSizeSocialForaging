data {
  int<lower=1> J;
  int<lower=1> Nt_t;
  int<lower=1> Nt_g;
  int<lower=1> Ng;

  int<lower=1> N_z1;
  int<lower=1> N_z2;

  // individual IDs
  array[N_z1] int<lower=1,upper=J> ind_z1;
  array[N_z2] int<lower=1,upper=J> ind_z2;

  // Assay IDs
  array[N_z1] int<lower=1,upper=Nt_t> Trial_t;
  array[N_z2] int<lower=1,upper=Nt_g> Trial_g;
  
  // group IDs
  array[N_z1] int<lower=1,upper=Ng> Group_t;
  array[N_z2] int<lower=1,upper=Ng> Group_g;

  // design variables
  vector[N_z1] Day;
  vector[N_z1] Seq_t;
  vector[N_z1] Grp_t;

  vector[N_z2] Seq_g;
  vector[N_z2] Grp_g;
  
  vector[N_z1] Events1;
  vector[N_z2] Events2;
  
  vector[N_z1] Events_opp1;
  vector[N_z2] Events_opp2;
  
  // focal counts
  array[N_z1] int<lower=0> Scr1;
  array[N_z2] int<lower=0> Scr2;

  array[N_z1] int<lower=0> n_events1;
  array[N_z2] int<lower=0> n_events2;

  // partner counts
  array[N_z1] int<lower=0> Scr_opp1;
  array[N_z2] int<lower=0> Scr_opp2;

  array[N_z1] int<lower=0> For_opp1;
  array[N_z2] int<lower=0> For_opp2;

  // response
  vector[N_z1] Rel_fit1;
  vector[N_z2] Rel_fit2;
}

transformed data {
  vector[N_z1] Prop_true1;
  vector[N_z2] Prop_true2;

  vector[N_z1] Prop_opp_true1;
  vector[N_z2] Prop_opp_true2;

  vector[N_z1] Prop_c1;
  vector[N_z2] Prop_c2;

  vector[N_z1] Prop_opp_c1;
  vector[N_z2] Prop_opp_c2;

  for(i in 1:N_z1){
    Prop_true1[i]     = Scr1[i] * 1.0 / n_events1[i];
    Prop_opp_true1[i] = Scr_opp1[i] * 1.0 / For_opp1[i];
  }

  for(i in 1:N_z2){
    Prop_true2[i]     = Scr2[i] * 1.0 / n_events2[i];
    Prop_opp_true2[i] = Scr_opp2[i] * 1.0 / For_opp2[i];
  }

  Prop_c1 = Prop_true1 - mean(Prop_true1);
  Prop_c2 = Prop_true2 - mean(Prop_true2);

  Prop_opp_c1 = Prop_opp_true1 - mean(Prop_opp_true1);
  Prop_opp_c2 = Prop_opp_true2 - mean(Prop_opp_true2);
}

parameters {

  //-----------------------------
  // fixed effects triad
  //-----------------------------

  real mu01;

  vector[7] beta1;

  real beta_main1_1;
  real beta_main2_1;

  real beta_quad1_1;
  real beta_quad2_1;

  real beta_inter1;

  real beta_events1;
  real beta_events2;
  
  real beta_events_opp1;
  real beta_events_opp2;
  //-----------------------------
  // fixed effects group
  //-----------------------------

  real mu02;

  vector[3] beta2;

  real beta_main1_2;
  real beta_main2_2;

  real beta_quad1_2;
  real beta_quad2_2;

  real beta_inter2;

  //-----------------------------
  // individual covariance matrix
  //-----------------------------

  matrix[J,2] zI;

  vector<lower=0>[2] sigma_I;
  cholesky_factor_corr[2] L_I;

  //-----------------------------
  // Assay effects
  //-----------------------------

  vector[Nt_t] zTrial1;
  vector[Nt_g] zTrial2;

  real<lower=0> sigma_Trial1;
  real<lower=0> sigma_Trial2;


  //-----------------------------
  // group effects
  //-----------------------------

  vector[Ng] zGroup1;
  vector[Ng] zGroup2;

  real<lower=0> sigma_Group1;
  real<lower=0> sigma_Group2;

  real<lower=0> sigma_e1;
  real<lower=0> sigma_e2;
}


transformed parameters{

  matrix[J,2] I;

  vector[Ng] gr1;
  vector[Ng] gr2;
  
  vector[Nt_t] tr1;
  vector[Nt_g] tr2;

  I = zI * diag_pre_multiply(sigma_I, L_I)';

  tr1 = sigma_Trial1 * zTrial1;
  tr2 = sigma_Trial2 * zTrial2;

  gr1 = sigma_Group1 * zGroup1;
  gr2 = sigma_Group2 * zGroup2;
}

model {

  vector[N_z1] mu1;
  vector[N_z2] mu2;

  //-----------------------------
  // triad predictor
  //-----------------------------

  mu1 =
      mu01
    + I[ind_z1,1]
    + beta1[1] * Day
    + beta1[2] * Seq_t
    + beta1[3] * Grp_t
    + beta1[4] * Day .* Seq_t
    + beta1[5] * Day .* Grp_t
    + beta1[6] * Seq_t .* Grp_t
    + beta1[7] * Day .* Seq_t .* Grp_t
    + beta_events1 * Events1 
    + beta_events_opp1 * Events_opp1
    + beta_main1_1 * Prop_c1
    + beta_main2_1 * Prop_opp_c1
    + beta_quad1_1 * (Prop_c1 .* Prop_c1)
    + beta_quad2_1 * (Prop_opp_c1 .* Prop_opp_c1)
    + beta_inter1 * (Prop_c1 .* Prop_opp_c1)
    + tr1[Trial_t]
    + gr1[Group_t];

  //-----------------------------
  // group predictor
  //-----------------------------

  mu2 =
      mu02
    + I[ind_z2,2]
    + beta2[1] * Seq_g
    + beta2[2] * Grp_g
    + beta2[3] * (Seq_g .* Grp_g)
    + beta_events2 * Events2
    + beta_events_opp2 * Events_opp2
    + beta_main1_2 * Prop_c2
    + beta_main2_2 * Prop_opp_c2
    + beta_quad1_2 * (Prop_c2 .* Prop_c2)
    + beta_quad2_2 * (Prop_opp_c2 .* Prop_opp_c2)
    + beta_inter2 * (Prop_c2 .* Prop_opp_c2)
    + tr2[Trial_g]
    + gr2[Group_g];

  //-----------------------------
  // likelihood
  //-----------------------------

  Rel_fit1 ~ normal(mu1, sigma_e1);
  Rel_fit2 ~ normal(mu2, sigma_e2);

  //-----------------------------
  // priors
  //-----------------------------

  to_vector(zI) ~ std_normal();
  L_I ~ lkj_corr_cholesky(2);

  sigma_I ~ exponential(2);

  zTrial1 ~ std_normal();
  zTrial2 ~ std_normal();

  sigma_Trial1 ~ exponential(2);
  sigma_Trial2 ~ exponential(2);
  
  zGroup1 ~ std_normal();
  zGroup2 ~ std_normal();

  sigma_Group1 ~ exponential(2);
  sigma_Group2 ~ exponential(2);

  sigma_e1 ~ exponential(1);
  sigma_e2 ~ exponential(1);

  mu01 ~ normal(0,1);
  mu02 ~ normal(0,1);

  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  
  beta_events1 ~ normal(0,1);
  beta_events2 ~ normal(0,1);
  
  beta_events_opp1 ~ normal(0,1);
  beta_events_opp2 ~ normal(0,1);
  
  beta_main1_1 ~ normal(0,1);
  beta_main2_1 ~ normal(0,1);
  beta_quad1_1 ~ normal(0,1);
  beta_quad2_1 ~ normal(0,1);
  beta_inter1 ~ normal(0,1);

  beta_main1_2 ~ normal(0,1);
  beta_main2_2 ~ normal(0,1);
  beta_quad1_2 ~ normal(0,1);
  beta_quad2_2 ~ normal(0,1);
  beta_inter2 ~ normal(0,1);
}

generated quantities {

  corr_matrix[2] Cor_I = multiply_lower_tri_self_transpose(L_I);

  real cor_int1_int2 = Cor_I[1,2]; // cross-context individual correlation

 // means on original scale
  real mean_feed1  = mu01;
  real mean_feed2  = mu02;
  
  real beta_focal1 = beta_main1_1;
  real beta_opp1   = beta_main2_1;

  real beta_focal2 = beta_main1_2;
  real beta_opp2   = beta_main2_2;

  real beta_int1   = beta_inter1;
  real beta_int2   = beta_inter2;

  vector[7] b1 = beta1;
  vector[3] b2 = beta2;

  real beta_focal_sq1 = beta_quad1_1*2;
  real beta_opp_sq1   = beta_quad2_1*2;

  real beta_focal_sq2 = beta_quad1_2*2;
  real beta_opp_sq2   = beta_quad2_2*2;

  // variances on original scale
  real V_I1 = square(sigma_I[1]);
  real V_I2 = square(sigma_I[2]);
  
  real V_Trial1 = square(sigma_Trial1);
  real V_Trial2 = square(sigma_Trial2);

  real V_Group1 = square(sigma_Group1);
  real V_Group2 = square(sigma_Group2);

  real V_e1 = square(sigma_e1);
  real V_e2 = square(sigma_e2);

  real V_total1 = V_I1 + V_Group1 + V_Trial1 + V_e1; // total variance triad
  real V_total2 = V_I2 + V_Group2 + V_Trial2 + V_e2; // total variance group

  real VC_I1 = V_I1 / V_total1;         // repeatability triad
  real VC_I2 = V_I2 / V_total2;         // repeatability group

  real VC_Trial1 = V_Trial1 / V_total1;         // repeatability triad
  real VC_Trial2 = V_Trial2 / V_total2;         // repeatability Trial
  
  real VC_Group1 = V_Group1 / V_total1;         // repeatability triad
  real VC_Group2 = V_Group2 / V_total2;         // repeatability group
  
  real VC_e1 = V_e1 / V_total1;         // repeatability triad
  real VC_e2 = V_e2 / V_total2;         // repeatability group

  real mean_prop1 = mean(Prop_true1);
  real mean_prop2 = mean(Prop_true2);

  real mean_opp1 = mean(Prop_opp_true1);
  real mean_opp2 = mean(Prop_opp_true2);

  real V_prop1 = variance(Prop_true1);
  real V_prop2 = variance(Prop_true2);

  real V_opp1 = variance(Prop_opp_true1);
  real V_opp2 = variance(Prop_opp_true2);
}


