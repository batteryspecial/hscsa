# =============================================================================
# Day 3: Statistical Hypothesis Testing
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("   HAMMETT SUBSTITUENT CONSTANT ANALYSIS - DAY 3\n")
cat("   Statistical Hypothesis Testing\n")
cat("============================================================\n")

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("\n[1] LOADING DATA...\n")

hammett <- read.csv("data/hammett_day_2_r.csv", stringsAsFactors = FALSE)

# Convert factors
hammett$heteroatom <- as.factor(hammett$heteroatom)
hammett$effect_type <- as.factor(hammett$effect_type)
hammett$functional_class <- as.factor(hammett$functional_class)

cat(sprintf("Loaded %d substituents\n", nrow(hammett)))

# -----------------------------------------------------------------------------
# 2. NORMALITY TESTING (Prerequisite for parametric tests)
# -----------------------------------------------------------------------------

cat("\n[2] NORMALITY TESTING...\n")
cat("H0: Data is normally distributed\n")
cat("H1: Data is not normally distributed\n")
cat("α = 0.05\n\n")

# Shapiro-Wilk test (best for small samples n < 50)
# Fortunately R does all the work for us
shapiro_meta <- shapiro.test(hammett$sigma_meta)
shapiro_para <- shapiro.test(hammett$sigma_para)
shapiro_res <- shapiro.test(hammett$sigma_resonance)

cat("--- Shapiro-Wilk Normality Tests ---\n")
cat(sprintf("σ_meta:      W = %.4f, p = %.4f %s\n", 
            shapiro_meta$statistic, shapiro_meta$p.value,
            ifelse(shapiro_meta$p.value < 0.05, "← NOT NORMAL", "← Normal")))
cat(sprintf("σ_para:      W = %.4f, p = %.4f %s\n", 
            shapiro_para$statistic, shapiro_para$p.value,
            ifelse(shapiro_para$p.value < 0.05, "← NOT NORMAL", "← Normal")))
cat(sprintf("σ_resonance: W = %.4f, p = %.4f %s\n", 
            shapiro_res$statistic, shapiro_res$p.value,
            ifelse(shapiro_res$p.value < 0.05, "← NOT NORMAL", "← Normal")))

cat("\nInterpretation: If p > 0.05, we fail to reject H0 (data is approximately normal).\n")
cat("Parametric tests (t-test, ANOVA) assume normality.\n")

# Visual normality check
png("output/day-3/11_normality_check.png", width = 1000, height = 400)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Q-Q plots
qqnorm(hammett$sigma_meta, main = expression("Q-Q Plot: " * sigma[meta]), 
       pch = 19, col = "#2166AC")
qqline(hammett$sigma_meta, col = "red", lwd = 2)

qqnorm(hammett$sigma_para, main = expression("Q-Q Plot: " * sigma[para]),
       pch = 19, col = "#2166AC")
qqline(hammett$sigma_para, col = "red", lwd = 2)

qqnorm(hammett$sigma_resonance, main = expression("Q-Q Plot: " * sigma[resonance]),
       pch = 19, col = "#2166AC")
qqline(hammett$sigma_resonance, col = "red", lwd = 2)

dev.off()
cat("\nSaved: output/day-3/11_normality_check.png\n")

# -----------------------------------------------------------------------------
# 3. PAIRED T-TEST: σ_meta vs σ_para
# -----------------------------------------------------------------------------

cat("\n[3] PAIRED T-TEST: σ_meta vs σ_para\n")
cat("============================================\n")
cat("Research Question: Are meta and para sigma values significantly different?\n\n")

cat("H0: μ(σ_meta) = μ(σ_para)  [no difference]\n")
cat("H1: μ(σ_meta) ≠ μ(σ_para)  [two-tailed]\n")
cat("α = 0.05\n\n")

# Perform paired t-test
# Paired because each substituent has BOTH a meta and para value
t_test_result <- t.test(hammett$sigma_meta, hammett$sigma_para, paired = TRUE)

cat("--- Results ---\n")
cat(sprintf("Mean σ_meta:        %.4f\n", mean(hammett$sigma_meta)))
cat(sprintf("Mean σ_para:        %.4f\n", mean(hammett$sigma_para)))
cat(sprintf("Mean difference:    %.4f\n", mean(hammett$sigma_meta - hammett$sigma_para)))
cat(sprintf("Standard Deviation: %.4f\n", sd(hammett$sigma_meta - hammett$sigma_para)))
cat(sprintf("t-statistic:        %.4f\n", t_test_result$statistic))
cat(sprintf("Degrees of freedom: %d\n", t_test_result$parameter))
cat(sprintf("p-value:            %.6f\n", t_test_result$p.value))
cat(sprintf("95%% CI:            [%.4f, %.4f]\n", 
            t_test_result$conf.int[1], t_test_result$conf.int[2]))

if (t_test_result$p.value < 0.05) {
  cat("\n→ CONCLUSION: Reject H0. σ_meta and σ_para are significantly different (p < 0.05).\n")
  cat("  This confirms that resonance effects create a systematic difference between positions.\n")
} else {
  cat("\n→ CONCLUSION: Fail to reject H0. No significant difference detected.\n")
}

# Effect size: Cohen's d for paired samples
diff <- hammett$sigma_meta - hammett$sigma_para
cohens_d <- mean(diff) / sd(diff)
cat(sprintf("\nEffect size (Cohen's d): %.3f ", cohens_d))
if (abs(cohens_d) < 0.2) {
  cat("[negligible]\n")
} else if (abs(cohens_d) < 0.5) {
  cat("[small]\n")
} else if (abs(cohens_d) < 0.8) {
  cat("[medium]\n")
} else {
  cat("[large]\n")
}

# Visualization
png("output/day-3/12_paired_ttest.png", width = 900, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Panel 1: Paired values connected by lines
plot(rep(1, nrow(hammett)), hammett$sigma_meta, 
     xlim = c(0.5, 2.5), ylim = range(c(hammett$sigma_meta, hammett$sigma_para)),
     pch = 19, col = "#2166AC", cex = 1.5,
     xlab = "", ylab = "σ value", xaxt = "n",
     main = "Paired Comparison: σ_meta vs σ_para")
points(rep(2, nrow(hammett)), hammett$sigma_para, pch = 19, col = "#B2182B", cex = 1.5)

# Connect pairs with lines
for (i in seq_len(nrow(hammett))) {
  line_col <- ifelse(hammett$sigma_para[i] < hammett$sigma_meta[i], "#2166AC50", "#B2182B50")
  lines(c(1, 2), c(hammett$sigma_meta[i], hammett$sigma_para[i]), col = line_col, lwd = 1.5)
}

axis(1, at = c(1, 2), labels = c(expression(sigma[meta]), expression(sigma[para])))
abline(h = 0, lty = 2, col = "gray50")

# Add means
points(1, mean(hammett$sigma_meta), pch = 18, col = "black", cex = 3)
points(2, mean(hammett$sigma_para), pch = 18, col = "black", cex = 3)
legend("topright", legend = c("Individual", "Mean"), pch = c(19, 18), 
       col = c("gray50", "black"), bty = "n")

# Panel 2: Distribution of differences
hist(diff, breaks = 12, col = "#67A9CF", border = "white",
     main = expression("Distribution of Differences (" * sigma[meta] - sigma[para] * ")"),
     xlab = expression(sigma[meta] - sigma[para]))
abline(v = 0, lty = 2, col = "red", lwd = 2)
abline(v = mean(diff), lty = 1, col = "blue", lwd = 2)
legend("topright", legend = c("Zero", paste("Mean =", round(mean(diff), 3))),
       col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n")

dev.off()
cat("Saved: output/day-3/12_paired_ttest.png\n")

# -----------------------------------------------------------------------------
# 4. ONE-WAY ANOVA: σ_para by Heteroatom Class
# -----------------------------------------------------------------------------

cat("\n[4] ONE-WAY ANOVA: σ_para by Heteroatom\n")
cat("============================================\n")
cat("Research Question: Do different heteroatoms have different mean σ_para values?\n\n")

cat("H0: μ_C = μ_N = μ_O = μ_halogen = ... [all means equal]\n")
cat("H1: At least one group mean differs\n")
cat("α = 0.05\n\n")

# Group heteroatoms for sufficient sample size
hammett$hetero_group <- as.character(hammett$heteroatom)
hammett$hetero_group[hammett$hetero_group %in% c("F", "Cl", "Br", "I")] <- "Halogen"
hammett$hetero_group[hammett$hetero_group %in% c("Si")] <- "Other"
hammett$hetero_group <- as.factor(hammett$hetero_group)

cat("--- Group Sizes ---\n")
print(table(hammett$hetero_group))

cat("\n--- Group Means ---\n")
group_means <- aggregate(sigma_para ~ hetero_group, data = hammett, 
                         FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
print(group_means)

# Perform ANOVA
anova_result <- aov(sigma_para ~ hetero_group, data = hammett)
anova_summary <- summary(anova_result)

cat("\n--- ANOVA Table ---\n")
print(anova_summary)

# Extract F-statistic and p-value
f_stat <- anova_summary[[1]]$`F value`[1]
p_val <- anova_summary[[1]]$`Pr(>F)`[1]

cat(sprintf("\nF-statistic: %.4f\n", f_stat))
cat(sprintf("p-value:     %.6f\n", p_val))

if (p_val < 0.05) {
  cat("\n→ CONCLUSION: Reject H0. Significant differences exist between heteroatom groups (p < 0.05).\n")
  
  # Post-hoc test: Tukey's HSD
  cat("\n--- Post-hoc Analysis: Tukey's HSD ---\n")
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
  
  # Identify significant pairs
  tukey_df <- as.data.frame(tukey_result$hetero_group)
  sig_pairs <- rownames(tukey_df)[tukey_df$`p adj` < 0.05]
  if (length(sig_pairs) > 0) {
    cat("\nSignificant pairwise differences (p < 0.05):\n")
    for (pair in sig_pairs) {
      cat(sprintf("  • %s\n", pair))
    }
  }
} else {
  cat("\n→ CONCLUSION: Fail to reject H0. No significant differences between groups.\n")
}

# ANOVA assumptions check
cat("\n--- Checking ANOVA Assumptions ---\n")

# Levene's test for homogeneity of variance (using Bartlett's test)
# Filter to groups with n >= 2
groups_with_data <- names(table(hammett$hetero_group))[table(hammett$hetero_group) >= 2]
hammett_filtered <- hammett[hammett$hetero_group %in% groups_with_data, ]

if (length(unique(hammett_filtered$hetero_group)) >= 2) {
  bartlett_result <- bartlett.test(sigma_para ~ hetero_group, data = hammett_filtered)
  cat(sprintf("Bartlett's test for homogeneity of variance: p = %.4f %s\n",
              bartlett_result$p.value,
              ifelse(bartlett_result$p.value < 0.05, "← Variances differ!", "← OK")))
} else {
  cat("Insufficient groups for Bartlett's test\n")
}

# Visualization
png("output/day-3/13_anova_heteroatom.png", width = 900, height = 600)
par(mar = c(6, 5, 4, 2))

# Boxplot with individual points
boxplot(sigma_para ~ hetero_group, data = hammett,
        col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
        ylab = expression(sigma[para]),
        main = expression("ANOVA: " * sigma[para] * " by Heteroatom Group"),
        cex.lab = 1.2, cex.main = 1.3)

# Add individual points with jitter
set.seed(42)
for (i in seq_along(levels(hammett$hetero_group))) {
  group <- levels(hammett$hetero_group)[i]
  y_vals <- hammett$sigma_para[hammett$hetero_group == group]
  x_vals <- jitter(rep(i, length(y_vals)), amount = 0.15)
  points(x_vals, y_vals, pch = 19, col = "#00000080", cex = 1.2)
}

abline(h = 0, lty = 2, col = "gray50")

# Add significance annotation if p < 0.05
if (p_val < 0.05) {
  mtext(sprintf("ANOVA: F = %.2f, p = %.4f *", f_stat, p_val), 
        side = 3, line = 0.3, cex = 0.9)
} else {
  mtext(sprintf("ANOVA: F = %.2f, p = %.4f (n.s.)", f_stat, p_val), 
        side = 3, line = 0.3, cex = 0.9)
}

dev.off()
cat("Saved: output/day-3/13_anova_heteroatom.png\n")

# -----------------------------------------------------------------------------
# 5. ONE-WAY ANOVA: σ_para by Effect Type
# -----------------------------------------------------------------------------

cat("\n[5] ONE-WAY ANOVA: σ_para by Effect Type\n")
cat("============================================\n")

anova_effect <- aov(sigma_para ~ effect_type, data = hammett)
anova_effect_summary <- summary(anova_effect)

cat("\n--- ANOVA Table ---\n")
print(anova_effect_summary)

f_stat2 <- anova_effect_summary[[1]]$`F value`[1]
p_val2 <- anova_effect_summary[[1]]$`Pr(>F)`[1]

cat(sprintf("\nF-statistic: %.4f\n", f_stat2))
cat(sprintf("p-value:     %.2e\n", p_val2))

if (p_val2 < 0.05) {
  cat("\n→ CONCLUSION: Reject H0. Effect types have significantly different σ_para values.\n")
  cat("  (This validates our EDG/EWG classification!)\n")
}

# Effect size: eta-squared
ss_between <- anova_effect_summary[[1]]$`Sum Sq`[1]
ss_total <- sum(anova_effect_summary[[1]]$`Sum Sq`)
eta_squared <- ss_between / ss_total
cat(sprintf("\nEffect size (η²): %.3f (%.1f%% of variance explained by effect type)\n", 
            eta_squared, eta_squared * 100))

cat("\n--- Tukey HSD Post-hoc Test ---\n")
tukey_effect <- TukeyHSD(anova_effect)
print(tukey_effect)

png("output/day-3/13_anova_effect_type.png", width = 1000, height = 600)

par(mar = c(7, 5, 4, 2))
effect_order <- aggregate(sigma_para ~ effect_type, data = hammett, FUN = mean)
effect_order <- effect_order[order(effect_order$sigma_para), ]
hammett$effect_type_ordered <- factor(hammett$effect_type, levels = effect_order$effect_type)

# Color palette: blue for EDG, red for EWG, gray for neutral
effect_colors <- c("EDG" = "#2166AC", "weak_EDG" = "#67A9CF", 
                   "neutral" = "#999999",
                   "weak_EWG" = "#EF8A62", "EWG" = "#B2182B")
box_colors <- effect_colors[levels(hammett$effect_type_ordered)]

boxplot(sigma_para ~ effect_type_ordered, data = hammett,
        col = box_colors, border = "black",
        ylab = expression(sigma[para]),
        xlab = "",
        main = expression("ANOVA: " * sigma[para] * " by Electronic Effect Type"),
        cex.lab = 1.2, cex.main = 1.3, las = 2)

# Add individual points with jitter
set.seed(42)
for (i in seq_along(levels(hammett$effect_type_ordered))) {
  group <- levels(hammett$effect_type_ordered)[i]
  y_vals <- hammett$sigma_para[hammett$effect_type_ordered == group]
  x_vals <- jitter(rep(i, length(y_vals)), amount = 0.15)
  points(x_vals, y_vals, pch = 19, col = "#00000080", cex = 1.2)
}

abline(h = 0, lty = 2, col = "gray50", lwd = 1.5)

# Add significance annotation
mtext(sprintf("F = %.2f, p = %.2e ***, η² = %.2f", f_stat2, p_val2, eta_squared), 
      side = 3, line = 0.3, cex = 0.95)

# Add group means as diamonds
group_means_effect <- aggregate(sigma_para ~ effect_type_ordered, data = hammett, FUN = mean)
points(seq_len(nrow(group_means_effect)), group_means_effect$sigma_para, 
       pch = 18, col = "black", cex = 2.5)

# Legend
legend("topleft", legend = c("EDG (electron-donating)", "Weak EDG", 
       "Neutral", "Weak EWG", "EWG (electron-withdrawing)", "Group mean"),
       fill = c(effect_colors, NA),
       pch = c(NA, NA, NA, NA, NA, 18),
       col = c(rep(NA, 5), "black"),
       border = c(rep("black", 5), NA),
       bty = "n", cex = 0.85)

dev.off()
cat("Saved: output/day-3/13_anova_effect_type.png\n")

# -----------------------------------------------------------------------------
# 6. LINEAR REGRESSION: Predicting σ_para from σ_meta
# -----------------------------------------------------------------------------

cat("\n[6] LINEAR REGRESSION: σ_para ~ σ_meta\n")
cat("============================================\n")
cat("Research Question: How well does σ_meta predict σ_para?\n\n")

# Fit linear model
lm_model <- lm(sigma_para ~ sigma_meta, data = hammett)
lm_summary <- summary(lm_model)

cat("--- Model: σ_para = β₀ + β₁·σ_meta + ε ---\n\n")
cat("Coefficients:\n")
print(lm_summary$coefficients)

cat(sprintf("\n--- Model Fit ---\n"))
cat(sprintf("R²:           %.4f (%.1f%% variance explained)\n", 
            lm_summary$r.squared, lm_summary$r.squared * 100))
cat(sprintf("Adjusted R²:  %.4f\n", lm_summary$adj.r.squared))
cat(sprintf("Residual SE:  %.4f\n", lm_summary$sigma))
cat(sprintf("F-statistic:  %.2f (p = %.2e)\n", 
            lm_summary$fstatistic[1], 
            pf(lm_summary$fstatistic[1], lm_summary$fstatistic[2], 
               lm_summary$fstatistic[3], lower.tail = FALSE)))

# Interpretation
cat("\n--- Interpretation ---\n")
cat(sprintf("Intercept (β₀ = %.3f): Predicted σ_para when σ_meta = 0\n", coef(lm_model)[1]))
cat(sprintf("Slope (β₁ = %.3f): For each unit increase in σ_meta, σ_para increases by %.3f\n",
            coef(lm_model)[2], coef(lm_model)[2]))
cat(sprintf("\nThe slope > 1 suggests para position amplifies electronic effects.\n"))

# Residual analysis
hammett$predicted <- predict(lm_model)
hammett$residuals <- residuals(lm_model)

cat("\n--- Residual Analysis ---\n")
cat("Large residuals indicate substituents where resonance dominates:\n\n")

# Top positive residuals (σ_para higher than predicted = resonance withdrawal)
cat("Positive residuals (σ_para > predicted, resonance acceptors):\n")
pos_resid <- hammett[order(-hammett$residuals), c("substituent", "sigma_meta", "sigma_para", "residuals")]
print(head(pos_resid, 5), row.names = FALSE)

# Top negative residuals (σ_para lower than predicted = resonance donation)
cat("\nNegative residuals (σ_para < predicted, resonance donors):\n")
neg_resid <- hammett[order(hammett$residuals), c("substituent", "sigma_meta", "sigma_para", "residuals")]
print(head(neg_resid, 5), row.names = FALSE)

# Visualization
png("output/day-3/14_regression_analysis.png", width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Panel 1: Regression plot
plot(hammett$sigma_meta, hammett$sigma_para,
     pch = 19, col = "#2166AC", cex = 1.5,
     xlab = expression(sigma[meta]),
     ylab = expression(sigma[para]),
     main = expression("Linear Regression: " * sigma[para] * " ~ " * sigma[meta]))

abline(lm_model, col = "#B2182B", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 1.5)  # y = x reference

# Label outliers (large residuals)
outlier_idx <- which(abs(hammett$residuals) > 0.25)
text(hammett$sigma_meta[outlier_idx], hammett$sigma_para[outlier_idx],
     labels = hammett$substituent[outlier_idx], pos = 4, cex = 0.7, col = "gray30")

legend("topleft", 
       legend = c(sprintf("Regression (R² = %.2f)", lm_summary$r.squared),
                  "y = x (no resonance)"),
       col = c("#B2182B", "gray50"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.9)

# Panel 2: Residuals vs Fitted
plot(hammett$predicted, hammett$residuals,
     pch = 19, col = "#2166AC", cex = 1.5,
     xlab = expression("Fitted " * sigma[para]),
     ylab = "Residuals",
     main = "Residuals vs Fitted Values")

abline(h = 0, lty = 2, col = "red", lwd = 2)

# Label outliers
text(hammett$predicted[outlier_idx], hammett$residuals[outlier_idx],
     labels = hammett$substituent[outlier_idx], pos = 4, cex = 0.7, col = "gray30")

# Add loess smoother to check for patterns
if (nrow(hammett) > 10) {
  loess_fit <- loess(residuals ~ predicted, data = hammett)
  pred_seq <- seq(min(hammett$predicted), max(hammett$predicted), length.out = 50)
  loess_pred <- predict(loess_fit, newdata = data.frame(predicted = pred_seq))
  lines(pred_seq, loess_pred, col = "#4DAF4A", lwd = 2)
  legend("topright", legend = "LOESS smoother", col = "#4DAF4A", lwd = 2, bty = "n")
}

dev.off()
cat("\nSaved: output/day-3/14_regression_analysis.png\n")

# Normality of residuals
shapiro_resid <- shapiro.test(hammett$residuals)
cat(sprintf("\nShapiro-Wilk test on residuals: W = %.4f, p = %.4f %s\n",
            shapiro_resid$statistic, shapiro_resid$p.value,
            ifelse(shapiro_resid$p.value > 0.05, "← Normal", "← Not normal")))

# -----------------------------------------------------------------------------
# 7. CORRELATION TESTS: Pearson vs Spearman
# -----------------------------------------------------------------------------

cat("\n[7] CORRELATION TESTS\n")
cat("============================================\n")
cat("Comparing Pearson (parametric) vs Spearman (non-parametric) correlations\n\n")

# Variables to test
vars_to_test <- c("sigma_I", "pi", "MR", "Es")

cat("--- Correlations with σ_para ---\n\n")
cat(sprintf("%-12s %10s %10s %10s %10s\n", "Variable", "Pearson r", "p-value", "Spearman ρ", "p-value"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (var in vars_to_test) {
  if (var %in% colnames(hammett)) {
    # Remove NAs for this variable
    complete_idx <- !is.na(hammett[[var]])
    
    # Pearson
    pearson <- cor.test(hammett$sigma_para[complete_idx], hammett[[var]][complete_idx], 
                        method = "pearson")
    # Spearman
    spearman <- cor.test(hammett$sigma_para[complete_idx], hammett[[var]][complete_idx], 
                         method = "spearman")
    
    cat(sprintf("%-12s %10.3f %10.4f %10.3f %10.4f\n", 
                var, pearson$estimate, pearson$p.value,
                spearman$estimate, spearman$p.value))
  }
}

cat("\nInterpretation:\n")
cat("• Pearson measures linear relationships (assumes normality)\n")
cat("• Spearman measures monotonic relationships (rank-based, robust to outliers)\n")
cat("• Similar values suggest linear relationship; different values suggest non-linearity\n")

# -----------------------------------------------------------------------------
# 8. MULTIPLE REGRESSION: σ_para ~ σ_meta + σ_I
# -----------------------------------------------------------------------------

cat("\n[8] MULTIPLE REGRESSION\n")
cat("============================================\n")
cat("Can we improve prediction by adding σ_I?\n\n")

# Remove rows with NA in sigma_I
hammett_complete <- hammett[!is.na(hammett$sigma_I), ]

# Simple model
model1 <- lm(sigma_para ~ sigma_meta, data = hammett_complete)

# Multiple regression
model2 <- lm(sigma_para ~ sigma_meta + sigma_I, data = hammett_complete)

cat("--- Model Comparison ---\n")
cat(sprintf("Model 1 (σ_meta only):       R² = %.4f, Adj.R² = %.4f\n",
            summary(model1)$r.squared, summary(model1)$adj.r.squared))
cat(sprintf("Model 2 (σ_meta + σ_I):      R² = %.4f, Adj.R² = %.4f\n",
            summary(model2)$r.squared, summary(model2)$adj.r.squared))

# F-test for model comparison
anova_compare <- anova(model1, model2)
cat("\n--- ANOVA: Model Comparison ---\n")
print(anova_compare)

if (anova_compare$`Pr(>F)`[2] < 0.05) {
  cat("\n→ Adding σ_I significantly improves the model.\n")
} else {
  cat("\n→ Adding σ_I does NOT significantly improve the model.\n")
  cat("  (σ_meta already captures the inductive information)\n")
}

# -----------------------------------------------------------------------------
# 9. SAVE RESULTS
# -----------------------------------------------------------------------------

cat("\n[9] SAVING RESULTS...\n")

# Save updated dataset
write.csv(hammett, "data/hammett_day_3_r.csv", row.names = FALSE)
cat("Saved: data/hammett_day_3_r.csv\n")

# Save summary
sink("docs/day-3/day3_summary.txt")
cat("============================================================\n")
cat("HAMMETT CONSTANT ANALYSIS - DAY 3 SUMMARY\n")
cat("Statistical Hypothesis Testing\n")
cat("============================================================\n\n")

cat("1. NORMALITY TESTS (Shapiro-Wilk)\n")
cat(sprintf("   σ_meta:      p = %.4f %s\n", shapiro_meta$p.value,
            ifelse(shapiro_meta$p.value > 0.05, "(normal)", "(not normal)")))
cat(sprintf("   σ_para:      p = %.4f %s\n", shapiro_para$p.value,
            ifelse(shapiro_para$p.value > 0.05, "(normal)", "(not normal)")))

cat("\n2. PAIRED T-TEST: σ_meta vs σ_para\n")
cat(sprintf("   t = %.4f, df = %d, p = %.6f\n", 
            t_test_result$statistic, t_test_result$parameter, t_test_result$p.value))
cat(sprintf("   Mean difference: %.4f\n", mean(diff)))
cat(sprintf("   Cohen's d: %.3f\n", cohens_d))
cat(sprintf("   Conclusion: %s\n", 
            ifelse(t_test_result$p.value < 0.05, "Significantly different", "Not significantly different")))

cat("\n3. ANOVA: σ_para by Heteroatom\n")
cat(sprintf("   F = %.4f, p = %.6f\n", f_stat, p_val))
cat(sprintf("   Conclusion: %s\n",
            ifelse(p_val < 0.05, "Groups differ significantly", "No significant difference")))

cat("\n4. ANOVA: σ_para by Effect Type\n")
cat(sprintf("   F = %.4f, p = %.2e\n", f_stat2, p_val2))
cat(sprintf("   η² = %.3f (%.1f%% variance explained)\n", eta_squared, eta_squared * 100))

cat("\n5. LINEAR REGRESSION: σ_para ~ σ_meta\n")
cat(sprintf("   R² = %.4f\n", lm_summary$r.squared))
cat(sprintf("   Equation: σ_para = %.3f + %.3f × σ_meta\n", coef(lm_model)[1], coef(lm_model)[2]))
cat(sprintf("   Residual SE: %.4f\n", lm_summary$sigma))

cat("\n6. KEY RESIDUALS (Resonance-dominated substituents)\n")
cat("   Positive (resonance acceptors): NO, COCH3, CF3\n")
cat("   Negative (resonance donors): N(CH3)2, NH2, OH\n")

sink()
cat("Saved: docs/day-3/day_3_summary.txt\n")

# -----------------------------------------------------------------------------
# COMPLETION
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("                  DAY 3 COMPLETE!\n")
cat("============================================================\n")
