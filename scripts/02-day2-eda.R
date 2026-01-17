# nolint file: line_length_linter

# =============================================================================
# Hammett Substituent Constant Analysis
# Day 2: Exploratory Data Analysis & Correlation Studies
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("   HAMMETT SUBSTITUENT CONSTANT ANALYSIS - DAY 2\n")
cat("   Exploratory Data Analysis & Correlation Studies\n")
cat("============================================================\n")

# -----------------------------------------------------------------------------
# 1. LOAD CLEANED DATA FROM DAY 1
# -----------------------------------------------------------------------------

cat("\n[1] LOADING DATA...\n")

hammett <- read.csv("data/hammett_clean.csv", stringsAsFactors = FALSE)
hammett_ext <- read.csv("data/hammett_extended.csv", stringsAsFactors = FALSE)

# Convert factors
hammett$heteroatom <- as.factor(hammett$heteroatom)
hammett$effect_type <- as.factor(hammett$effect_type)
hammett$functional_class <- as.factor(hammett$functional_class)

cat(sprintf("Loaded %d substituents with %d variables\n", nrow(hammett), ncol(hammett)))

# -----------------------------------------------------------------------------
# 2. CORRELATION MATRIX ANALYSIS
# -----------------------------------------------------------------------------

cat("\n[2] CORRELATION ANALYSIS...\n")

# Select numeric columns for correlation
numeric_cols <- c("sigma_meta", "sigma_para", "sigma_I", "sigma_v", "pi", "Es", "MR") 
numeric_data <- hammett[, numeric_cols]

# Calculate correlation matrix (handling NAs)
cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

cat("\n--- Correlation Matrix ---\n")
print(round(cor_matrix, 3))

# Key correlations to highlight
cat("\n--- Key Correlations ---\n")
cat(sprintf("σ_meta vs σ_para: r = %.3f (inductive+resonance relationship)\n",  cor_matrix["sigma_meta", "sigma_para"])) # nolint: line_length_linter.
cat(sprintf("σ_meta vs σ_I: r = %.3f (meta ≈ inductive effect)\n", cor_matrix["sigma_meta", "sigma_I"]))
cat(sprintf("σ_para vs σ_I: r = %.3f (para includes resonance beyond inductive)\n", cor_matrix["sigma_para", "sigma_I"]))# nolint: line_length_linter.
cat(sprintf("σ_para vs π (hydrophobicity): r = %.3f\n", cor_matrix["sigma_para", "pi"]))
cat(sprintf("MR vs Es (size parameters):   r = %.3f\n", cor_matrix["MR", "Es"]))

# -----------------------------------------------------------------------------
# 3. CORRELATION HEATMAP VISUALIZATION
# -----------------------------------------------------------------------------

cat("\n[3] GENERATING CORRELATION HEATMAP...\n")
png("output/day-2/06_correlation_heatmap.png", width = 900, height = 700)

# Use layout to properly separate heatmap and legend
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))

# panel 1 is heatmap
par(mar = c(6, 6, 4, 1))

# Create color palette (blue = negative, white = zero, red = positive)
n_colors <- 100
color_palette <- colorRampPalette(c("#2166AC", "#67A9CF",
                                    "#D1E5F0", "white",
                                    "#FDDBC7", "#EF8A62",
                                    "#B2182B"))(n_colors)

# Plot heatmap
image(seq_len(ncol(cor_matrix)), seq_len(nrow(cor_matrix)), t(cor_matrix)
      [, rev(seq_len(nrow(cor_matrix)))],
      col = color_palette, zlim = c(-1, 1),
      xlab = "", ylab = "", axes = FALSE,
      main = "Correlation Matrix: Hammett & Molecular Descriptors")

# Add axis labels
axis(1, at = seq_len(ncol(cor_matrix)), labels = colnames(cor_matrix), las = 2, cex.axis = 0.9)
axis(2, at = seq_len(nrow(cor_matrix)), labels = rev(rownames(cor_matrix)), las = 1, cex.axis = 0.9)

# Add correlation values as text
for (i in seq_len(ncol(cor_matrix))) {
  for (j in seq_len(nrow(cor_matrix))) {
    val <- cor_matrix[nrow(cor_matrix) - j + 1, i]
    if (!is.na(val)) {
      text_col <- ifelse(abs(val) > 0.5, "white", "black")
      text(i, j, sprintf("%.2f", val), cex = 0.7, col = text_col)
    }
  }
}

# Panel 2 is color legend
par(mar = c(6, 0, 4, 3))
legend_seq <- seq(-1, 1, length.out = n_colors)
image(1, legend_seq, matrix(legend_seq, nrow = 1), col = color_palette,
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(4, at = c(-1, -0.5, 0, 0.5, 1), las = 1, cex.axis = 0.8)
par(fig = c(0.85, 0.95, 0.15, 0.85), new = TRUE, mar = c(0, 0, 0, 0))
mtext("r", side = 4, line = 2.5, cex = 0.9)

dev.off()
cat("Saved: output/day-2/06_correlation_heatmap.png\n")

# -----------------------------------------------------------------------------
# 4. RESONANCE vs INDUCTIVE DECOMPOSITION
# -----------------------------------------------------------------------------

cat("\n[4] RESONANCE vs INDUCTIVE ANALYSIS...\n")

# Calculate resonance contribution
# The difference between sigma_para and sigma_meta reveals resonance contribution
# σ_para = σ_inductive + σ_resonance
# σ_meta ≈ σ_inductive (resonance blocked at meta)
# Therefore: σ_resonance ≈ σ_para - σ_meta

hammett$sigma_resonance <- hammett$sigma_para - hammett$sigma_meta

cat("\n--- Resonance Contribution (σ_para - σ_meta) ---\n")
cat(sprintf("Mean resonance effect: %.3f\n", mean(hammett$sigma_resonance)))
cat(sprintf("SD resonance effect:   %.3f\n", sd(hammett$sigma_resonance)))
cat(sprintf("Range: [%.3f, %.3f]\n", min(hammett$sigma_resonance),
            max(hammett$sigma_resonance)))

# Classify resonance behavior
hammett$resonance_type <- ifelse(hammett$sigma_resonance < -0.15,
                                 "Strong +R (donor)",
                          ifelse(hammett$sigma_resonance > 0.15,
                                 "Strong -R (acceptor)",
                                 "Weak/No resonance"))

cat("\n--- Resonance Classification ---\n")
print(table(hammett$resonance_type))

# Top resonance donors (most negative sigma_resonance)
cat("\n--- Top Resonance Donors (+R groups) ---\n")
donors <- hammett[order(hammett$sigma_resonance), 
                  c("substituent", "sigma_meta", "sigma_para", "sigma_resonance")]
print(head(donors, 8), row.names = FALSE)

# Top resonance acceptors (most positive sigma_resonance)
cat("\n--- Top Resonance Acceptors (-R groups) ---\n")
acceptors <- hammett[order(-hammett$sigma_resonance),
                     c("substituent", "sigma_meta", "sigma_para", "sigma_resonance")]
print(head(acceptors, 5), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. RESONANCE CONTRIBUTION VISUALIZATION
# -----------------------------------------------------------------------------

cat("\n[5] GENERATING RESONANCE ANALYSIS PLOT...\n")

png("output/day-2/07_resonance_contribution.png", width = 900, height = 700)
par(mar = c(8, 5, 4, 2))

# Sort by resonance contribution
hammett_sorted <- hammett[order(hammett$sigma_resonance), ]

# Create bar colors based on resonance type
bar_colors <- ifelse(hammett_sorted$sigma_resonance < -0.15, "#2166AC",  # donors = blue
               ifelse(hammett_sorted$sigma_resonance > 0.15, "#B2182B",   # acceptors = red
                      "#999999"))                                          # neutral = gray

barplot(hammett_sorted$sigma_resonance, 
        names.arg = hammett_sorted$substituent,
        col = bar_colors, border = NA,
        ylab = expression("Resonance Contribution (" * sigma[para] - sigma[meta] * ")"),
        main = "Resonance Effect by Substituent",
        las = 2, cex.names = 0.8, cex.lab = 1.2)

abline(h = 0, lty = 1, col = "black", lwd = 1)
abline(h = c(-0.15, 0.15), lty = 2, col = "gray50", lwd = 1)

legend("topleft", 
       legend = c("Strong +R (electron donor)", "Weak/None", "Strong -R (electron acceptor)"),
       fill = c("#2166AC", "#999999", "#B2182B"),
       bty = "n", cex = 0.9)

mtext("Negative = resonance donation | Positive = resonance withdrawal", 
      side = 3, line = 0.3, cex = 0.9)

dev.off()
cat("Saved: output/day-2/07_resonance_contribution.png\n")

# -----------------------------------------------------------------------------
# 6. FUNCTIONAL GROUP CLASS ANALYSIS
# -----------------------------------------------------------------------------

cat("\n[6] FUNCTIONAL GROUP CLASS ANALYSIS...\n")

# Summary statistics by functional class
func_stats <- aggregate(cbind(sigma_meta, sigma_para, sigma_resonance) ~ functional_class, 
                        data = hammett, 
                        FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))

# Flatten the result
func_summary <- data.frame(
  functional_class = func_stats$functional_class,
  n = func_stats$sigma_meta[, "n"],
  mean_sigma_m = round(func_stats$sigma_meta[, "mean"], 3),
  mean_sigma_p = round(func_stats$sigma_para[, "mean"], 3),
  mean_resonance = round(func_stats$sigma_resonance[, "mean"], 3)
)
func_summary <- func_summary[order(func_summary$mean_sigma_p), ]

cat("\n--- Statistics by Functional Class ---\n")
print(func_summary, row.names = FALSE)

# Visualization
png("output/day-2/08_functional_class_comparison.png", width = 1000, height = 600)
par(mar = c(8, 5, 4, 8), xpd = TRUE)

# Prepare data
classes <- func_summary$functional_class
n_classes <- length(classes)
x_positions <- barplot(rep(0, n_classes), names.arg = classes, las = 2, ylim = c(-0.9, 0.9), plot = FALSE)

# Plot sigma_meta and sigma_para as grouped bars
bar_width <- 0.35
barplot(func_summary$mean_sigma_m, 
        col = "#67A9CF", border = NA,
        ylim = c(-0.9, 0.9),
        names.arg = classes, las = 2,
        ylab = expression("Mean " * sigma * " value"),
        main = "Electronic Effects by Functional Class",
        cex.names = 0.85, cex.lab = 1.2)

# Overlay sigma_para as points
points(x_positions, func_summary$mean_sigma_p, 
       pch = 18, col = "#B2182B", cex = 2)

abline(h = 0, lty = 2, col = "gray50", lwd = 1)

legend("topright", inset = c(-0.15, 0),
       legend = c(expression(sigma[meta] ~ "(bars)"),
                  expression(sigma[para] ~ "(diamonds)")),
       fill = c("#67A9CF", NA),
       pch = c(NA, 18),
       col = c(NA, "#B2182B"),
       pt.cex = 1.5, bty = "n", cex = 0.9)

dev.off()
cat("Saved: output/day-2/08_functional_class_comparison.png\n")

# -----------------------------------------------------------------------------
# 7. OUTLIER DETECTION (IQR METHOD)
# -----------------------------------------------------------------------------

cat("\n[7] OUTLIER DETECTION...\n")

# Function to detect outliers using IQR method
detect_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr

  return (x < lower | x > upper)
}

# Detect outliers for each sigma type
hammett$outlier_sigma_m <- detect_outliers(hammett$sigma_meta)
hammett$outlier_sigma_p <- detect_outliers(hammett$sigma_para)
hammett$outlier_resonance <- detect_outliers(hammett$sigma_resonance)

cat("\n--- Outliers in σ_meta ---\n")
outliers_m <- hammett[hammett$outlier_sigma_m,
                      c("substituent", "sigma_meta", "effect_type")]

if (nrow(outliers_m) > 0)
  print(outliers_m, row.names = FALSE) else cat("None\n")

cat("\n--- Outliers in σ_para ---\n")
outliers_p <- hammett[hammett$outlier_sigma_p, 
                      c("substituent", "sigma_para", "effect_type")]

if (nrow(outliers_p) > 0)
  print(outliers_p, row.names = FALSE) else cat("None\n")

cat("\n--- Outliers in Resonance Contribution ---\n")
outliers_r <- hammett[hammett$outlier_resonance, 
                      c("substituent", "sigma_resonance", "resonance_type")]

if (nrow(outliers_r) > 0) 
  print(outliers_r, row.names = FALSE) else cat("None\n")

# -----------------------------------------------------------------------------
# 8. EXTENDED DATASET ANALYSIS (σ+ and σ-)
# -----------------------------------------------------------------------------

cat("\n[8] EXTENDED ANALYSIS: σ+ and σ- CONSTANTS...\n")

# Clean extended dataset
hammett_ext$sigma_plus <- as.numeric(hammett_ext$sigma_plus)
hammett_ext$sigma_minus <- as.numeric(hammett_ext$sigma_minus)
hammett_ext_clean <- hammett_ext[!is.na(hammett_ext$sigma_plus) & !is.na(hammett_ext$sigma_minus),]

cat(sprintf("\nExtended dataset: %d substituents with σ+ and σ- values\n", 
            nrow(hammett_ext_clean)))

# σ+ is used for reactions with positive charge buildup (electrophilic)
# σ- is used for reactions with negative charge buildup (nucleophilic)

cat("\n--- σ+ Constants (Cation-stabilizing ability) ---\n")
cat("More negative σ+ = better at stabilizing positive charge\n")
cat("\nBest cation stabilizers:\n")
cation_stab <- hammett_ext_clean[order(hammett_ext_clean$sigma_plus), 
                                  c("substituent", "sigma_plus", "sigma_para")]
print(head(cation_stab, 8), row.names = FALSE)

cat("\n--- σ- Constants (Anion-stabilizing ability) ---\n")
cat("More positive σ- = better at stabilizing negative charge\n")
cat("\nBest anion stabilizers:\n")
anion_stab <- hammett_ext_clean[order(-hammett_ext_clean$sigma_minus),
                                 c("substituent", "sigma_minus", "sigma_para")]
print(head(anion_stab, 8), row.names = FALSE)

# Visualization: σ+ vs σ- comparison
png("output/day-2/09_sigma_plus_minus.png", width = 900, height = 700)
par(mar = c(5, 5, 4, 2))

plot(hammett_ext_clean$sigma_plus, hammett_ext_clean$sigma_minus,
     pch = 19, cex = 1.5,
     col = ifelse(hammett_ext_clean$effect_type == "EDG", "#2166AC",
           ifelse(hammett_ext_clean$effect_type == "EWG", "#B2182B", "#666666")),
     xlab = expression(sigma^"+"),
     ylab = expression(sigma^"-"),
     main = expression("Hammett " * sigma^"+" * " vs " * sigma^"-" * " Constants"),
     cex.lab = 1.3, cex.main = 1.4)

# Add diagonal line
abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)

# Label key points
key_subs <- c("NH2", "N(CH3)2", "OCH3", "NO2", "CN", "F", "Cl")
for (sub in key_subs) {
  idx <- which(hammett_ext_clean$substituent == sub)
  if (length(idx) > 0) {
    text(hammett_ext_clean$sigma_plus[idx], hammett_ext_clean$sigma_minus[idx],
         labels = sub, pos = 4, cex = 0.8, col = "gray30")
  }
}

# Add reference lines
abline(h = 0, v = 0, lty = 3, col = "gray70")

legend("topleft",
       legend = c("EDG", "EWG", "Other", expression(sigma^"+" == sigma^"-")),
       col = c("#2166AC", "#B2182B", "#666666", "gray50"),
       pch = c(19, 19, 19, NA), lty = c(NA, NA, NA, 2),
       bty = "n", cex = 0.9)

mtext(expression("Points below line: " * sigma^"-" > sigma^"+" * " (enhanced EWG for anions)"), 
      side = 3, line = 0.3, cex = 0.9)

dev.off()
cat("Saved: output/day-2/09_sigma_plus_minus.png\n")

# -----------------------------------------------------------------------------
# 9. PAIRWISE SCATTER PLOT MATRIX
# -----------------------------------------------------------------------------

cat("\n[9] GENERATING SCATTER PLOT MATRIX...\n")

png("output/day-2/10_scatter_matrix.png", width = 900, height = 900)

# Select key variables
plot_vars <- hammett[, c("sigma_meta", "sigma_para", "sigma_I", "sigma_resonance")]
colnames(plot_vars) <- c("σ_meta", "σ_para", "σ_I", "σ_resonance")

# Custom panel function for scatter plots
panel_scatter <- function(x, y, ...) {
  points(x, y, pch = 19, col = "#2166AC80", cex = 1.2)
  abline(lm(y ~ x), col = "#B2182B", lwd = 2)
}

# Custom panel function for correlation values
panel_cor <- function(x, y, ...) {
  # usr <- par("usr")
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs")
  txt <- sprintf("r = %.2f", r)
  text(0.5, 0.5, txt, cex = 1.5, col = ifelse(abs(r) > 0.7, "#B2182B", "black"))
}

pairs(plot_vars,
      lower.panel = panel_scatter,
      upper.panel = panel_cor,
      diag.panel = NULL,
      main = "Pairwise Relationships: Hammett Constants",
      cex.main = 1.3)

dev.off()
cat("Saved: output/day-2/10_scatter_matrix.png\n")

# -----------------------------------------------------------------------------
# 10. SAVE UPDATED DATA AND SUMMARY
# -----------------------------------------------------------------------------

cat("\n[10] SAVING RESULTS...\n")

# Save updated dataset with resonance calculations
write.csv(hammett, "data/hammett_day2.csv", row.names = FALSE)
cat("Saved: data/hammett_day2.csv\n")

# Save Day 2 summary
sink("docs/day-2-summary.txt")
cat("============================================================\n")
cat("HAMMETT CONSTANT ANALYSIS - DAY 2 SUMMARY\n")
cat("============================================================\n\n")

cat("KEY CORRELATIONS:\n")
cat(sprintf("  σ_meta vs σ_para: r = %.3f\n", cor_matrix["sigma_meta", "sigma_para"]))
cat(sprintf("  σ_meta vs σ_I:    r = %.3f\n", cor_matrix["sigma_meta", "sigma_I"]))
cat(sprintf("  σ_para vs σ_I:    r = %.3f\n", cor_matrix["sigma_para", "sigma_I"]))

cat("\nRESONANCE CONTRIBUTION (σ_para - σ_meta):\n")
cat(sprintf("  Mean: %.3f\n", mean(hammett$sigma_resonance)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(hammett$sigma_resonance), max(hammett$sigma_resonance)))

cat("\nRESONANCE CLASSIFICATION:\n")
print(table(hammett$resonance_type))

cat("\nTOP RESONANCE DONORS:\n")
print(head(donors, 5), row.names = FALSE)

cat("\nTOP RESONANCE ACCEPTORS:\n")
print(head(acceptors, 5), row.names = FALSE)

cat("\nFUNCTIONAL CLASS SUMMARY:\n")
print(func_summary, row.names = FALSE)

cat("\nOUTLIERS DETECTED:\n")
cat(sprintf("  In σ_meta: %d\n", sum(hammett$outlier_sigma_m)))
cat(sprintf("  In σ_para: %d\n", sum(hammett$outlier_sigma_p)))
cat(sprintf("  In resonance: %d\n", sum(hammett$outlier_resonance)))
sink()
cat("Saved: output/day-2/day2_summary.txt\n")

# -----------------------------------------------------------------------------
# COMPLETION
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("                  DAY 2 COMPLETE!\n")
cat("============================================================\n")
cat("\nGenerated Files:\n")
cat("  DATA:\n")
cat("    - data/hammett_day2.csv (with resonance calculations)\n")
cat("  PLOTS:\n")
cat("    - output/day-2/06_correlation_heatmap.png\n")
cat("    - output/day-2/07_resonance_contribution.png\n")
cat("    - output/day-2/08_functional_class_comparison.png\n")
cat("    - output/day-2/09_sigma_plus_minus.png\n")
cat("    - output/day-2/10_scatter_matrix.png\n")
cat("  SUMMARY:\n")
cat("    - docs/day2_summary.txt\n")
cat("\nKEY INSIGHTS:\n")
cat("  - σ_meta correlates strongly with σ_I (inductive effect)\n")
cat("  - Resonance contribution = σ_para - σ_meta\n")
cat("  - Amino groups are strongest resonance donors\n")
cat("  - Carbonyl/nitro groups are resonance acceptors\n")
cat("\nNEXT (Day 3): Statistical hypothesis testing\n")
cat("  - Paired t-test: σ_meta vs σ_para\n")
cat("  - ANOVA: Effect types and heteroatom classes\n")
cat("  - Normality tests\n")
cat("============================================================\n")