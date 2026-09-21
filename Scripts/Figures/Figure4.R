# Run Scripts/Analyses/Analysis_mod1.R and Scripts/Analyses/Analysis_mod2.R 
# to obtain the required model output results

library(plot3D)

md <- readRDS("Model_output/Mod2_results.RDS")
md1 <- readRDS("Model_output/Mod1_results.RDS")
post_mat <- as.data.frame(as.matrix(md))

#--------------------------------------------------
# Posterior mean coefficients
#--------------------------------------------------

# Triadic assays
b0_1  <- mean(post_mat$mu01)
b1_1  <- mean(post_mat$beta_main1_1)
b2_1  <- mean(post_mat$beta_main2_1)
b11_1 <- mean(post_mat$beta_quad1_1)
b22_1 <- mean(post_mat$beta_quad2_1)
b12_1 <- mean(post_mat$beta_inter1)

mean_prop1 <- mean(post_mat$mean_prop1)
mean_opp1  <- mean(post_mat$mean_opp1)


# Group assays
b0_2  <- mean(post_mat$mu02)
b1_2  <- mean(post_mat$beta_main1_2)
b2_2  <- mean(post_mat$beta_main2_2)
b11_2 <- mean(post_mat$beta_quad1_2)
b22_2 <- mean(post_mat$beta_quad2_2)
b12_2 <- mean(post_mat$beta_inter2)

mean_prop2 <- mean(post_mat$mean_prop2)
mean_opp2  <- mean(post_mat$mean_opp2)


#--------------------------------------------------
# Prediction grid
#--------------------------------------------------

partner <- seq(0, 1, length.out = 25)
focal   <- seq(0, 1, length.out = 25)

grid <- expand.grid(
  Partner = partner,
  Prop    = focal
)


#--------------------------------------------------
# Triadic surface
#--------------------------------------------------

Prop_c1    <- grid$Prop    - mean_prop1
Partner_c1 <- grid$Partner - mean_opp1


Fitness1 <-
  b0_1 +
  b1_1  * Prop_c1 +
  b2_1  * Partner_c1 +
  b11_1 * Prop_c1^2 +
  b22_1 * Partner_c1^2 +
  b12_1 * Prop_c1 * Partner_c1

zmat1 <- matrix(
  Fitness1,
  nrow = length(partner),
  ncol = length(focal),
  byrow = FALSE
)


#--------------------------------------------------
# Group surface
#--------------------------------------------------

Prop_c2    <- grid$Prop - mean_prop2
Partner_c2 <- grid$Partner - mean_opp2

Fitness2 <-
  b0_2 +
  b1_2  * Prop_c2 +
  b2_2  * Partner_c2 +
  b11_2 * Prop_c2^2 +
  b22_2 * Partner_c2^2 +
  b12_2 * Prop_c2 * Partner_c2

zmat2 <- matrix(
  Fitness2,
  nrow = length(partner),
  ncol = length(focal),
  byrow = FALSE
)

#--------------------------------------------------
# Common z-axis and black-and-white palette
#--------------------------------------------------

zlim <- range(c(zmat1, zmat2))

bw_palette <- colorRampPalette(
  c("black", "white")
)(25)


#--------------------------------------------------
# Plot
#--------------------------------------------------

par(
  mfrow = c(1, 2),
  mar = c(3, 3, 2, 2)
)


# Triadic assays
persp3D(
  x = partner,
  y = focal,
  z = zmat1,
  col = bw_palette,
  border = "black",
  alpha = 0.8,
  ticktype = "detailed",
  xlab = "",
  ylab = "",
  zlab="",
  theta = 30,
  phi = 50,
  xlim = c(0, 1),
  ylim = c(0, 1),
  zlim = c(0.4, 1.2),
  main = "Triadic assays",
  colkey = FALSE,
  cex.axis = 0.8
)

partner_br <- seq(0,1,length.out=200)

focal_br <-
  mean_prop1 -
  (
    b1_1 +
      b12_1*(partner_br - mean_opp1)
  )/(2*b11_1)

focal_br <- pmin(pmax(focal_br,0),1)

# Centred values
focal_c   <- focal_br   - mean_prop1
partner_c <- partner_br - mean_opp1

# Fitness at focal optimum
z_br <-
  b0_1 +
  b1_1  * focal_c +
  b2_1  * partner_c +
  b11_1 * focal_c^2 +
  b22_1 * partner_c^2 +
  b12_1 * focal_c * partner_c


lines3D(
  x = partner_br,
  y = focal_br,
  z = z_br + 0.015,
  add = TRUE,
  col = "red",
  lwd = 4
)

# Predicted fitness continued
z_br <-
  b0_1 + 
  b1_1  *(0 - mean_prop1) +
  b2_1  * (seq(max(partner_br),1,length.out=20) - mean_opp1) +
  b11_1 *  (0 - mean_prop1)^2 +
  b22_1 * (seq(max(partner_br),1,length.out=20) -mean_opp1)^2 +
  b12_1 * (0 - mean_prop1) * (seq(max(partner_br),1,length.out=20)-mean_opp1)

lines3D(
  x = seq(max(partner_br),1,length.out=20),
  y = rep(0.001,20),
  z = z_br+0.015,
  add = TRUE,
  col = "black",
  lwd = 4
)

eq <- seq(0, 1, length.out = 200)

eq_prop_c <- eq - mean_prop1
eq_opp_c  <- eq - mean_opp1

z_eq <-
  b0_1 +
  b1_1  * eq_prop_c +
  b2_1  * eq_opp_c +
  b11_1 * eq_prop_c^2 +
  b22_1 * eq_opp_c^2 +
  b12_1 * eq_prop_c * eq_opp_c

lines3D(
  x = eq,
  y = eq,
  z = z_eq+0.001,
  add = TRUE,
  col = "black",
  lwd = 2
)

#md1 <- readRDS("Mod3_noEIV_results.RDS")
### 1. CORRELATION PLOT WITH 95% CIs (BLUPs) ----

# Stan summary for Int_RE[J,2]
sm_RNj <- as.data.frame(rstan::summary(md1, pars = "Int_RE")$summary)
sm_RNj$par <- rownames(sm_RNj)

sm <- as.data.frame(rstan::summary(md1, pars = c("Int_RE", "mu0_1", "mu0_2", "beta_p1", "beta_p2"))$summary)
sm$par <- rownames(sm)


# population intercepts & slopes (means)
mu0_1   <- sm$mean[sm$par == "mu0_1"]
mu0_2   <- sm$mean[sm$par == "mu0_2"]
beta_p1 <- sm$mean[sm$par == "beta_p1"]
beta_p2 <- sm$mean[sm$par == "beta_p2"]

## Population reaction norm
x_pop <- seq(0, 1, length.out = 200)

y_pop <- plogis(
  mu0_1 +
    beta_p1 * (x_pop - mean_opp1)
)


## Height of the fitness surface
partner_pop_c <- x_pop - mean_opp1
focal_pop_c   <- y_pop - mean_prop1

z_pop <-
  b0_1 +
  b1_1  * focal_pop_c +
  b2_1  * partner_pop_c +
  b11_1 * focal_pop_c^2 +
  b22_1 * partner_pop_c^2 +
  b12_1 * focal_pop_c * partner_pop_c

## Draw the line
lines3D(
  x = x_pop,
  y = y_pop,
  z = z_pop + 0.015,
  add = TRUE,
  col = "black",
  lwd = 4
)



#--------------------------------------------------
# Group plot
#--------------------------------------------------
persp3D(
  x = partner,
  y = focal,
  z = zmat2,
  col = bw_palette,
  border = "black",
  alpha = 0.8,
  ticktype = "detailed",
  xlab = "",
  ylab = "",
  zlab="",
  theta = 30,
  phi = 50,
  xlim = c(0, 1),
  ylim = c(0, 1),
  zlim = c(0.4, 1.2),
  main = "Group assays",
  colkey = FALSE,
  cex.axis = 0.8
)

partner_br <- seq(0,1,length.out=200)

focal_br <-
  mean_prop2 -
  (
    b1_2 +
      b12_2*(partner_br - mean_opp2)
  )/(2*b11_2)

focal_br <- pmin(pmax(focal_br,0),1)

# Centred values
focal_c   <- focal_br   - mean_prop2
partner_c <- partner_br - mean_opp2

# Fitness at focal optimum
z_br <-
  b0_2 +
  b1_2  * focal_c +
  b2_2  * partner_c +
  b11_2 * focal_c^2 +
  b22_2 * partner_c^2 +
  b12_2 * focal_c * partner_c


eq <- seq(0, 1, length.out = 200)

eq_prop_c <- eq - mean_prop2
eq_opp_c  <- eq - mean_opp2

z_eq <-
  b0_2 +
  b1_2  * eq_prop_c +
  b2_2  * eq_opp_c +
  b11_2 * eq_prop_c^2 +
  b22_2 * eq_opp_c^2 +
  b12_2 * eq_prop_c * eq_opp_c

lines3D(
  x = eq,
  y = eq,
  z = z_eq + 0.001,
  add = TRUE,
  col = "black",
  lwd = 2
)


## Population reaction norm
x_pop <- seq(0, 1, length.out = 200)

y_pop <- plogis(
  mu0_2 +
    beta_p2 * (x_pop - mean_opp2)
)

partner_pop_c <- x_pop - mean_opp2
focal_pop_c   <- y_pop - mean_prop2

z_pop <-
  b0_2 +
  b1_2  * focal_pop_c +
  b2_2  * partner_pop_c +
  b11_2 * focal_pop_c^2 +
  b22_2 * partner_pop_c^2 +
  b12_2 * focal_pop_c * partner_pop_c

## Draw the line
lines3D(
  x = x_pop,
  y = y_pop,
  z = z_pop + 0.015,
  add = TRUE,
  col = "black",
  lwd = 4
)


ess1 <-
  (post_mat$mean_prop1 -
     post_mat$beta_main1_1 / (2 * post_mat$beta_quad1_1) +
     (post_mat$beta_inter1 * post_mat$mean_opp1) / (2 * post_mat$beta_quad1_1)) /
  (1 + post_mat$beta_inter1 / (2 * post_mat$beta_quad1_1) )

quantile(
  ess1,
  c(0.025, 0.5, 0.975)
)

ess2 <-
  (post_mat$mean_prop2 -
      post_mat$beta_main1_2 / (2 * post_mat$beta_quad1_2) +
      (post_mat$beta_inter2 * post_mat$mean_opp2) / (2 * post_mat$beta_quad1_2)) /  (1 +
      post_mat$beta_inter2 / (2 * post_mat$beta_quad1_2))

quantile(
  ess2,
  c(0.025, 0.5, 0.975)
)



