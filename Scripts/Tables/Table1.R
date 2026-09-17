library(ggplot2)

# ============================================================
# Load posterior samples
# ============================================================

md <- readRDS("Model_output/Mod2_results.RDS")

post <- as.data.frame(as.matrix(md))


# ============================================================
# Calculate adaptive-dynamics quantities for each posterior draw
# ============================================================

calc_stability <- function(post,
                           assay,
                           b1_name,
                           b11_name,
                           b12_name,
                           mean_prop_name,
                           mean_opp_name) {
  
  b1  <- post[[b1_name]]
  b11 <- post[[b11_name]]
  b12 <- post[[b12_name]]
  mp  <- post[[mean_prop_name]]
  mo  <- post[[mean_opp_name]]
  
  # Singular strategy
  s_star <- (
    2 * b11 * mp +
      b12 * mo -
      b1
  ) / (2 * b11 + b12)
  
  # Convergence stability:
  # 2*b11 + b12 < 0
  convergence <- 2 * b11 + b12
  
  # Evolutionary stability:
  # 2*b11 < 0
  evolutionary_stability <- 2 * b11
  
  data.frame(
    assay = assay,
    s_star = s_star,
    convergence = convergence,
    evolutionary_stability = evolutionary_stability
  )
}


# ============================================================
# Calculate for both assays
# ============================================================

triadic_stability <- calc_stability(
  post = post,
  assay = "Small group",
  b1_name = "beta_main1_1",
  b11_name = "beta_quad1_1",
  b12_name = "beta_inter1",
  mean_prop_name = "mean_prop1",
  mean_opp_name = "mean_opp1"
)

group_stability <- calc_stability(
  post = post,
  assay = "Large group",
  b1_name = "beta_main1_2",
  b11_name = "beta_quad1_2",
  b12_name = "beta_inter2",
  mean_prop_name = "mean_prop2",
  mean_opp_name = "mean_opp2"
)


# Combine posterior draws
stability <- rbind(
  triadic_stability,
  group_stability
)


# ============================================================
# Summarise posterior
# ============================================================

summarise_posterior <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  c(
    median = unname(q[2]),
    lower  = unname(q[1]),
    upper  = unname(q[3])
  )
}


stability_summary <- do.call(
  rbind,
  lapply(
    split(stability, stability$assay),
    function(x) {
      
      s  <- summarise_posterior(x$s_star)
      cs <- summarise_posterior(x$convergence)
      es <- summarise_posterior(x$evolutionary_stability)
      
      data.frame(
        assay = unique(x$assay),
        
        s_median  = s[1],
        s_lower   = s[2],
        s_upper   = s[3],
        
        cs_median = cs[1],
        cs_lower  = cs[2],
        cs_upper  = cs[3],
        
        es_median = es[1],
        es_lower  = es[2],
        es_upper  = es[3]
      )
    }
  )
)

rownames(stability_summary) <- NULL

stability_summary

summ_ci <- function(x) {
  paste0(
    sprintf("%.3f", x["median"]),
    " [",
    sprintf("%.3f", x["lower"]),
    ", ",
    sprintf("%.3f", x["upper"]),
    "]"
  )
}

```r
# ------------------------------------------------------------
# Format posterior median and 95% credible interval
# ------------------------------------------------------------

fmt_ci <- function(median, lower, upper) 
  { sprintf( "%.3f [%.3f, %.3f]", median, lower, upper ) }


# ------------------------------------------------------------
# LaTeX table: Adaptive dynamics
# ------------------------------------------------------------

latex_lines <- c(
  "\\begin{table}[ht]",
  "\\caption{Posterior estimates of the singular strategy and its evolutionary properties.}",
  "\\centering",
  "\\small",
  "\\begin{tabular}{lccc}",
  "\\hline",
  "\\textbf{Assay} & \\textbf{$s^*$} & \\textbf{Convergence stability} & \\textbf{Evolutionary stability} \\\\",
  "\\hline",
  
  paste0(
    stability_summary$assay[2], " & ",
    fmt_ci(
      stability_summary$s_median[2],
      stability_summary$s_lower[2],
      stability_summary$s_upper[2]
    ), " & ",
    fmt_ci(
      stability_summary$cs_median[2],
      stability_summary$cs_lower[2],
      stability_summary$cs_upper[2]
    ), " & ",
    fmt_ci(
      stability_summary$es_median[2],
      stability_summary$es_lower[2],
      stability_summary$es_upper[2]
    ),
    " \\\\"
  ),
  
  paste0(
    stability_summary$assay[1], " & ",
    fmt_ci(
      stability_summary$s_median[1],
      stability_summary$s_lower[1],
      stability_summary$s_upper[1]
    ), " & ",
    fmt_ci(
      stability_summary$cs_median[1],
      stability_summary$cs_lower[1],
      stability_summary$cs_upper[1]
    ), " & ",
    fmt_ci(
      stability_summary$es_median[1],
      stability_summary$es_lower[1],
      stability_summary$es_upper[1]
    ),
    " \\\\"
  ),
  
  "\\hline",
  "\\end{tabular}",
  "\\label{tab:adaptive_dynamics}",
  "\\end{table}"
)

cat(latex_lines, sep = "\n")
