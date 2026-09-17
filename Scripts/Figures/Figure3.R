library(dplyr)
library(tidyr)
library(ggplot2)

md <- readRDS("Model_output/Mod1_results.RDS")

### 1. CORRELATION PLOT WITH 95% CIs (BLUPs) ----

# Stan summary for Int_RE[J,2]
sm_RNj <- as.data.frame(rstan::summary(md, pars = "Int_RE")$summary)
sm_RNj$par <- rownames(sm_RNj)

sm <- as.data.frame(
  rstan::summary(
    md,
    pars = c("Int_RE", "mu0_1", "mu0_2")
  )$summary
)
sm$par <- rownames(sm)


# Population intercepts: POSTERIOR MEDIANS
mu0_1   <- sm$`50%`[sm$par == "mu0_1"]
mu0_2   <- sm$`50%`[sm$par == "mu0_2"]

# Extract ID and dimension from names like "Int_RE[12,1]"
sm_RNj$par <- rownames(sm_RNj)

blups <- sm_RNj %>%
  tidyr::extract(
    col = par,
    into = c("ID", "dim"),
    regex = "Int_RE\\[(\\d+),(\\d+)\\]",
    convert = TRUE
  ) %>%
  filter(dim %in% c(1, 2)) %>%
  mutate(
    context = if_else(dim == 1, "triadic", "group"),
    
    # POSTERIOR MEDIAN individual intercept
    intercept = if_else(
      dim == 1,
      mu0_1 + `50%`,
      mu0_2 + `50%`
    ),
    
    # 95% posterior intervals
    lower = if_else(
      dim == 1,
      mu0_1 + `2.5%`,
      mu0_2 + `2.5%`
    ),
    
    upper = if_else(
      dim == 1,
      mu0_1 + `97.5%`,
      mu0_2 + `97.5%`
    ),
    
    m = plogis(intercept),
    l = plogis(lower),
    u = plogis(upper)
  ) %>%
  select(ID, context, m, l, u) %>%
  pivot_wider(
    id_cols = ID,
    names_from = context,
    values_from = c(m, l, u),
    names_sep = "_"
  )


gg_cor_int <- ggplot(blups, aes(x = m_triadic, y = m_group)) +
  geom_point(alpha = 1) +
  geom_errorbar(
    aes(ymin = l_group, ymax = u_group),
    alpha = 0.2
  ) +
  geom_errorbarh(
    aes(xmin = l_triadic, xmax = u_triadic),
    alpha = 0.2
  ) +
  labs(
    title = "Individual foraging strategy across contexts",
    x = "Scrounging tendency \n(small group)",
    y = "Scrounging tendency \n(large group)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text  = element_text(colour = "black", size = 10),
    axis.title = element_text(colour = "black", size = 12)
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0.01, 0.01)
  )


gg_cor_int

