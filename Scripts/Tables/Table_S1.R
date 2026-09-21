# Run Scripts/Analyses/Analysis_mod1.R 
# to obtain the required model output results

md <- readRDS("Model_output/Mod1_results.RDS")

post_mat <- as.data.frame(as.matrix(md))

summ_ci <- function(x, digits = 3) {
  q <- unname(quantile(as.numeric(x), c(0.025, 0.5, 0.975), na.rm = TRUE))
  sprintf("%.*f\\;[%.*f,\\;%.*f]",
          digits, q[2], digits, q[1], digits, q[3])
}

diff_ci <- function(x, y, digits = 3) {
  q <- unname(quantile(as.numeric(x) - as.numeric(y),
                       c(0.025, 0.5, 0.975), na.rm = TRUE))
  sprintf("%.*f\\;[%.*f,\\;%.*f]",
          digits, q[2], digits, q[1], digits, q[3])
}

dash <- "--"

latex_lines <- c(
  "\\begin{table}[ht]",
  "\\caption{Posterior medians and 95\\% credible intervals}",
  "\\centering",
  "\\small",
  "\\begin{tabular}{lccc}",
  "\\hline",
  "\\textbf{Parameter} & \\textbf{Small group} & \\textbf{Large group} & \\textbf{Difference (S--L)} \\\\",
  "\\hline",
  
  "\\multicolumn{4}{l}{\\textbf{Fixed effects}} \\\\",
  paste0("$\\mu_0$ Intercept & ",
         summ_ci(post_mat$mu0_1), " & ",
         summ_ci(post_mat$mu0_2), " & ",
         diff_ci(post_mat$mu0_1, post_mat$mu0_2), " \\\\"),
  "\\multicolumn{4}{l}{\\,\\,\\,\\textit{Controls}} \\\\",
  paste0("$\\beta_1$ Test day (1st or 2nd) (Day) & ",
         summ_ci(post_mat[["beta_t[1]"]]), " & ", dash, " & ", dash, " \\\\"),
  paste0("$\\beta_2$ Assay order (Seq) & ",
         summ_ci(post_mat[["beta_t[2]"]]), " & ",
         summ_ci(post_mat[["beta_g[1]"]]), " & ",
         diff_ci(post_mat[["beta_t[2]"]], post_mat[["beta_g[1]"]]), " \\\\"),
  paste0("$\\beta_3$ Large group assay first (Grp) & ",
         summ_ci(post_mat[["beta_t[3]"]]), " & ",
         summ_ci(post_mat[["beta_g[2]"]]), " & ",
         diff_ci(post_mat[["beta_t[3]"]], post_mat[["beta_g[2]"]]), " \\\\"),
  paste0("$\\beta_4$ Day*Seq & ",
         summ_ci(post_mat[["beta_t[4]"]]), " & ", dash, " & ", dash, " \\\\"),
  paste0("$\\beta_5$ Day*Grp & ",
         summ_ci(post_mat[["beta_t[5]"]]), " & ", dash, " & ", dash, " \\\\"),
  paste0("$\\beta_6$ Seq*Grp & ",
         summ_ci(post_mat[["beta_t[6]"]]), " & ",
         summ_ci(post_mat[["beta_g[3]"]]), " & ",
         diff_ci(post_mat[["beta_t[6]"]], post_mat[["beta_g[3]"]]), " \\\\"),
  paste0("$\\beta_7$ Day*Seq*Grp & ",
         summ_ci(post_mat[["beta_t[7]"]]), " & ", dash, " & ", dash, " \\\\"),
  
  "\\multicolumn{4}{l}{\\,\\,\\,\\textit{Average response to partner tactic}} \\\\",
  paste0("$\\beta_8$ Partner scrounging & ",
         summ_ci(post_mat$beta_p1), " & ",
         summ_ci(post_mat$beta_p2), " & ",
         diff_ci(post_mat$beta_p1, post_mat$beta_p2), " \\\\"),
  "\\hline",
  "\\multicolumn{4}{l}{\\textbf{Random effect}} \\\\",
  paste0("$V_{\\text{ID}}$ & ",
         summ_ci(post_mat$V_int1), " & ",
         summ_ci(post_mat$V_int2), " & ",
         diff_ci(post_mat$V_int1, post_mat$V_int2), " \\\\"),
  paste0("$V_{\\text{AssayID}}$ & ",
         summ_ci(post_mat$V_trial1), " & ",
         summ_ci(post_mat$V_trial2), " & ",
         diff_ci(post_mat$V_trial1, post_mat$V_trial2), " \\\\"),
  paste0("$V_{\\text{GroupID}}$ & ",
         summ_ci(post_mat$V_group1), " & ",
         summ_ci(post_mat$V_group2), " & ",
         diff_ci(post_mat$V_group1, post_mat$V_group2), " \\\\"),
  paste0("$V_{\\text{Total}}$ & ",
         summ_ci(post_mat$V_total1), " & ",
         summ_ci(post_mat$V_total2), " & ",
         diff_ci(post_mat$V_total1, post_mat$V_total2), " \\\\"),
  "\\hline",
  "\\multicolumn{4}{l}{\\textbf{Variance components}} \\\\",
  paste0("$VC_{\\text{Among-individual}}$ & ",
         summ_ci(post_mat$VC_int1), " & ",
         summ_ci(post_mat$VC_int2), " & ",
         diff_ci(post_mat$VC_int1, post_mat$VC_int2), " \\\\"),
  paste0("$VC_{\\text{Social response}}$ & ",
         summ_ci(post_mat$VC_social1), " & ",
         summ_ci(post_mat$VC_social2), " & ",
         diff_ci(post_mat$VC_social1, post_mat$VC_social2), " \\\\"),
  paste0("$VC_{\\text{Residual}}$ & ",
         summ_ci(post_mat$VC_residual1), " & ",
         summ_ci(post_mat$VC_residual2), " & ",
         diff_ci(post_mat$VC_residual1, post_mat$VC_residual2), " \\\\"),
  paste0("$VC_{\\text{Within-individual}}$ & ",
         summ_ci(post_mat$VC_within1), " & ",
         summ_ci(post_mat$VC_within2), " & ",
         diff_ci(post_mat$VC_within1, post_mat$VC_within2), " \\\\"),
  paste0("$VC_{\\text{AssayID}}$ & ",
         summ_ci(post_mat$VC_trial1), " & ",
         summ_ci(post_mat$VC_trial2), " & ",
         diff_ci(post_mat$VC_trial1, post_mat$VC_trial2), " \\\\"),
  paste0("$VC_{\\text{GroupID}}$ & ",
         summ_ci(post_mat$VC_group1), " & ",
         summ_ci(post_mat$VC_group2), " & ",
         diff_ci(post_mat$VC_group1, post_mat$VC_group2), " \\\\"),
  "\\hline",
  "\\multicolumn{4}{l}{\\textbf{Correlations}} \\\\",
  paste0("$\\mathrm{cor}(\\text{ID}_{small},\\text{ID}_{large})$ & \\multicolumn{2}{c}{",
         summ_ci(post_mat$cor_int1_int2),
         "} & \\multicolumn{1}{c}{--} \\\\"),
  "\\hline",
  "\\end{tabular}",
  "\\label{tab:stan_results_full}",
  "\\end{table}"
)

cat(paste(latex_lines, collapse = "\n"))
