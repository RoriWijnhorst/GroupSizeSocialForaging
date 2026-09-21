# Run Scripts/Analyses/Analysis_mod1.R and Scripts/Analyses/Analysis_mod2.R 
# to obtain model output results

library(ggplot2)

# ============================================================
# Load posterior samples
# ============================================================

md <- readRDS("Model_output/Mod2_results.RDS")
md1 <- readRDS("Model_output/Mod1_results.RDS")

post1 <- as.data.frame(as.matrix(md))
post2 <- as.data.frame(as.matrix(md1))


# ============================================================
# Function to calculate and plot PIP
# ============================================================

make_pip <- function(
    post_mat,
    assay,
    b0_name,
    b1_name,
    b2_name,
    b11_name,
    b22_name,
    b12_name,
    mean_prop_name,
    mean_opp_name
) {
  
  # ----------------------------------------------------------
  # Posterior mean coefficients
  # ----------------------------------------------------------
  
  b0  <- mean(post_mat[[b0_name]])
  b1  <- mean(post_mat[[b1_name]])
  b2  <- mean(post_mat[[b2_name]])
  b11 <- mean(post_mat[[b11_name]])
  b22 <- mean(post_mat[[b22_name]])
  b12 <- mean(post_mat[[b12_name]])
  
  mean_prop <- mean(post_mat[[mean_prop_name]])
  mean_opp  <- mean(post_mat[[mean_opp_name]])
  
  
  # ----------------------------------------------------------
  # Grid of resident and alternative strategies
  # ----------------------------------------------------------
  
  pip <- expand.grid(
    s_R = seq(0, 1, length.out = 500),
    s_M = seq(0, 1, length.out = 500)
  )
  
  
  # ----------------------------------------------------------
  # Centre focal strategies
  # ----------------------------------------------------------
  
  pip$s_R_c <- pip$s_R - mean_prop
  pip$s_M_c <- pip$s_M - mean_prop
  
  
  # ----------------------------------------------------------
  # Alternative focal fitness
  # ----------------------------------------------------------
  
  pip$W_M <-
    b0 +
    b1  * pip$s_M_c +
    b2  * (pip$s_R - mean_opp) +
    b11 * pip$s_M_c^2 +
    b22 * (pip$s_R - mean_opp)^2 +
    b12 * pip$s_M_c * (pip$s_R - mean_opp)
  
  
  # ----------------------------------------------------------
  # Resident focal fitness
  # ----------------------------------------------------------
  
  pip$W_R <-
    b0 +
    b1  * pip$s_R_c +
    b2  * (pip$s_R - mean_opp) +
    b11 * pip$s_R_c^2 +
    b22 * (pip$s_R - mean_opp)^2 +
    b12 * pip$s_R_c * (pip$s_R - mean_opp)
  
  
  # ----------------------------------------------------------
  # Invasion fitness
  # ----------------------------------------------------------
  
  pip$invasion_fitness <- pip$W_M - pip$W_R
  
  pip$can_invade <- pip$invasion_fitness > 0
  
  
  # ----------------------------------------------------------
  # Singular / self-consistent strategy
  # ----------------------------------------------------------
  
  s_star <- (
    2 * b11 * mean_prop +
      b12 * mean_opp -
      b1
  ) / (
    2 * b11 + b12
  )
  
  
  # ----------------------------------------------------------
  # Best-response strategy
  # ----------------------------------------------------------
  
  pip$best_response <-
    mean_prop -
    (
      b1 +
        b12 * (pip$s_R - mean_opp)
    ) / (2 * b11)
  
  # Restrict to possible strategies
  pip$best_response <- pmin(
    pmax(pip$best_response, 0),
    1
  )
  
  
  # ----------------------------------------------------------
  # Check diagonal
  # ----------------------------------------------------------
  
  diag_check <- subset(
    pip,
    abs(s_M - s_R) < 0.003
  )
  
  cat("\n", assay, "\n")
  cat("Singular strategy:", s_star, "\n")
  cat(
    "Range of invasion fitness near diagonal:",
    range(diag_check$invasion_fitness),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    pip,
    aes(
      x = s_R,
      y = s_M
    )
  ) +
    
    # PIP
    geom_raster(
      aes(fill = can_invade)
    ) +
    
    # Resident = alternative
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "white",
      linewidth = 1.2
    ) +
    
    # Best-response function
    geom_line(
      aes(
        x = s_R,
        y = best_response
      ),
      color = "red",
      linewidth = 1.6
    ) +
    
    scale_fill_manual(
      values = c(
        "FALSE" = "grey5",
        "TRUE" = "grey80"
      ),
      labels = c(
        "FALSE" = "Focal does worse",
        "TRUE" = "Focal does better"
      )
    ) +
    
    scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    
    labs(
      title = assay,
      x = "Partner scrounging proportion",
      y = "Focal scrounging proportion",
      fill = NULL
    ) +
    
    coord_fixed(ratio = 1) +
    
    theme_classic(base_size = 13) +
    
    theme(
      legend.position = "right",
      legend.key.height = unit(0.8, "cm"),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
  
  
  # Return everything
  list(
    pip = pip,
    s_star = s_star,
    plot = p,
    coefficients = c(
      b0 = b0,
      b1 = b1,
      b2 = b2,
      b11 = b11,
      b22 = b22,
      b12 = b12,
      mean_prop = mean_prop,
      mean_opp = mean_opp
    )
  )
}


# ============================================================
# Run for both assays
# ============================================================

small <- make_pip(
  post_mat = post1,
  assay = "Small group",
  b0_name = "mu01",
  b1_name = "beta_main1_1",
  b2_name = "beta_main2_1",
  b11_name = "beta_quad1_1",
  b22_name = "beta_quad2_1",
  b12_name = "beta_inter1",
  mean_prop_name = "mean_prop1",
  mean_opp_name = "mean_opp1"
)


large <- make_pip(
  post_mat = post1,
  assay = "Large group",
  b0_name = "mu02",
  b1_name = "beta_main1_2",
  b2_name = "beta_main2_2",
  b11_name = "beta_quad1_2",
  b22_name = "beta_quad2_2",
  b12_name = "beta_inter2",
  mean_prop_name = "mean_prop2",
  mean_opp_name = "mean_opp2"
)


# ============================================================
# Display plots together
# ============================================================

library(patchwork)

combined_plot <- small$plot + large$plot +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

combined_plot



