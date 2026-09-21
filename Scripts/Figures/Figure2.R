# Run Scripts/Analyses/Analysis_mod1.R 
# to obtain model output results

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

md <- readRDS("Model_output/Mod1_results.RDS")

# Stan summary for Int_RE[J,2]
sm_RNj <- as.data.frame(rstan::summary(md, pars = "Int_RE")$summary)
sm_RNj$par <- rownames(sm_RNj)

sm <- as.data.frame(
  rstan::summary(
    md,
    pars = c("Int_RE", "mu0_1", "mu0_2", "beta_p1", "beta_p2")
  )$summary
)
sm$par <- rownames(sm)


# Population intercepts & slopes: POSTERIOR MEDIANS
mu0_1   <- sm$`50%`[sm$par == "mu0_1"]
mu0_2   <- sm$`50%`[sm$par == "mu0_2"]
beta_p1 <- sm$`50%`[sm$par == "beta_p1"]
beta_p2 <- sm$`50%`[sm$par == "beta_p2"]

### 2. REACTION NORMS ----

# Individual deviations from Int_RE[J,4]
# Use POSTERIOR MEDIANS
rn_medians <- sm %>%
  filter(grepl("^Int_RE\\[", par)) %>%
  tidyr::extract(
    par,
    into = c("ID", "dim"),
    regex = "Int_RE\\[(\\d+),(\\d+)\\]",
    convert = TRUE
  ) %>%
  select(ID, dim, `50%`) %>%
  tidyr::pivot_wider(
    names_from = dim,
    values_from = `50%`,
    names_prefix = "d"
  ) %>%
  transmute(
    ID,
    
    # logit-scale intercepts = population median + individual deviation
    int_t = mu0_1 + d1,
    int_g = mu0_2 + d2
  )



x_seq <- seq(0, 1, length.out = 100)  # raw partner proportion scale


# Social environment mean (mean_centering)
# Social environment mean
mean_triad <- mean(rstan::extract(md, pars = "x_obs1")$x_obs1[1, ])
mean_group <- mean(rstan::extract(md, pars = "x_obs2")$x_obs2[1, ])


# ----- TRIADIC LINES (reaction norms) -----
triadic_lines <- rn_medians %>%
  tidyr::expand_grid(x = x_seq) %>%
  mutate(
    x_c = x - mean_triad,
    eta = int_t + beta_p1 * x_c,
    p   = plogis(eta)
  )


# ----- GROUP LINES (reaction norms) -----
group_lines <- rn_medians %>%
  tidyr::expand_grid(x = x_seq) %>%
  mutate(
    x_c = x - mean_group,
    eta = int_g + beta_p2 * x_c,
    p   = plogis(eta)
  )


# ----- POPULATION LINES -----

pop_triadic <- data.frame(
  x = x_seq,
  p = plogis(
    mu0_1 + beta_p1 * (x_seq - mean_triad)
  )
)

pop_group <- data.frame(
  x = x_seq,
  p = plogis(
    mu0_2 + beta_p2 * (x_seq - mean_group)
  )
)


# ===== PLOT 1: TRIADIC REACTION NORMS =====

gg_triadic <- ggplot() +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "darkgrey"
  ) +
  geom_hline(
    yintercept = plogis(mu0_1),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_line(
    data = triadic_lines,
    aes(x = x, y = p, group = ID),
    color = "darkgrey",
    alpha = 0.4
  ) +
  geom_line(
    data = pop_triadic,
    aes(x = x, y = p),
    linewidth = 1.4,
    colour = "black",
    inherit.aes = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text  = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 12)
  ) +
  xlab("Partner scrounging proportion") +
  ylab("Focal scrounging proportion") +
  ggtitle("Small group reaction norms") +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  )


# ===== PLOT 2: GROUP REACTION NORMS =====

gg_group <- ggplot() +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "darkgrey"
  ) +
  geom_hline(
    yintercept = plogis(mu0_2),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_line(
    data = group_lines,
    aes(x = x, y = p, group = ID),
    color = "darkgrey",
    alpha = 0.4
  ) +
  geom_line(
    data = pop_group,
    aes(x = x, y = p),
    linewidth = 1.3,
    colour = "black",
    inherit.aes = FALSE
  ) +
  theme_classic() +
  theme(
    axis.text  = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 12)
  ) +
  xlab("Partner scrounging proportion") +
  ylab("Focal scrounging proportion") +
  ggtitle("Large group reaction norms") +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  )


# Layout: two reaction-norm panels 
(gg_triadic + gg_group)

