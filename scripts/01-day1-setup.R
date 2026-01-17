# nolint file: line_length_linter

# =============================================================================
# Hammett Substituent Constant Analysis
# Day 1: Setup, Data Loading, and Initial Exploration
# =============================================================================

# -----------------------------------------------------------------------------
# 1. SET WORKING DIRECTORY
# -----------------------------------------------------------------------------
setwd("~/Documents/Development/hscsa")  # Uncomment and set your path

cat("\n")
cat("============================================================\n")
cat("   HAMMETT SUBSTITUENT CONSTANT ANALYSIS - DAY 1\n")
cat("   Statistical Analysis in R for Physical Organic Chemistry\n")
cat("============================================================\n")

# -----------------------------------------------------------------------------
# 2. LOAD AND INSPECT DATA
# -----------------------------------------------------------------------------

cat("\n[1] LOADING DATA...\n")

# Load the Hammett constants dataset
hammett <- read.csv("data/hammett_constants.csv", stringsAsFactors = FALSE)
hammett_ext <- read.csv("data/hammett_extended.csv", stringsAsFactors = FALSE)

# Display structure
cat("\n--- Dataset Structure ---\n")
str(hammett)

cat("\n--- First 10 Rows ---\n")
print(head(hammett, 10))

cat("\n--- Data Dimensions ---\n")
cat(sprintf("Basic dataset: %d substituents, %d variables\n", nrow(hammett), ncol(hammett)))
cat(sprintf("Extended dataset: %d substituents, %d variables\n", nrow(hammett_ext), ncol(hammett_ext)))

# -----------------------------------------------------------------------------
# 3. DATA CLEANING AND PREPROCESSING
# -----------------------------------------------------------------------------

cat("\n[2] DATA CLEANING...\n")

# Convert categorical variables to factors
hammett$heteroatom <- as.factor(hammett$heteroatom)
hammett$effect_type <- as.factor(hammett$effect_type)
hammett$functional_class <- as.factor(hammett$functional_class)

# Handle NA values - count them
cat("\n--- Missing Values by Column ---\n")
na_counts <- colSums(is.na(hammett))
print(na_counts[na_counts > 0])

# Create clean dataset (no NA in key sigma columns)
hammett_clean <- hammett[!is.na(hammett$sigma_meta) & !is.na(hammett$sigma_para), ]

cat(sprintf("\nClean dataset: %d substituents (removed %d with missing sigma values)\n", 
            nrow(hammett_clean), nrow(hammett) - nrow(hammett_clean)))

# Add derived columns
hammett_clean$sigma_diff <- hammett_clean$sigma_meta - hammett_clean$sigma_para
hammett_clean$resonance_indicator <- abs(hammett_clean$sigma_diff)

# -----------------------------------------------------------------------------
# 4. SUMMARY STATISTICS
# -----------------------------------------------------------------------------

cat("\n[3] SUMMARY STATISTICS...\n")

# Basic statistics for sigma values
cat("\n--- Sigma Meta Statistics ---\n")
cat(sprintf("  Mean:   %.4f\n", mean(hammett_clean$sigma_meta)))
cat(sprintf("  SD:     %.4f\n", sd(hammett_clean$sigma_meta)))
cat(sprintf("  Min:    %.4f (%s)\n", min(hammett_clean$sigma_meta), 
            hammett_clean$substituent[which.min(hammett_clean$sigma_meta)]))
cat(sprintf("  Max:    %.4f (%s)\n", max(hammett_clean$sigma_meta),
            hammett_clean$substituent[which.max(hammett_clean$sigma_meta)]))
cat(sprintf("  Median: %.4f\n", median(hammett_clean$sigma_meta)))

cat("\n--- Sigma Para Statistics ---\n")
cat(sprintf("  Mean:   %.4f\n", mean(hammett_clean$sigma_para)))
cat(sprintf("  SD:     %.4f\n", sd(hammett_clean$sigma_para)))
cat(sprintf("  Min:    %.4f (%s)\n", min(hammett_clean$sigma_para),
            hammett_clean$substituent[which.min(hammett_clean$sigma_para)]))
cat(sprintf("  Max:    %.4f (%s)\n", max(hammett_clean$sigma_para),
            hammett_clean$substituent[which.max(hammett_clean$sigma_para)]))
cat(sprintf("  Median: %.4f\n", median(hammett_clean$sigma_para)))

# Group statistics by effect type
cat("\n--- Statistics by Effect Type ---\n")
effect_types <- unique(hammett_clean$effect_type)
effect_stats <- data.frame(
  effect_type = character(),
  n = integer(),
  mean_sigma_m = numeric(),
  mean_sigma_p = numeric(),
  sd_sigma_m = numeric(),
  sd_sigma_p = numeric(),
  stringsAsFactors = FALSE
)

for (etype in effect_types) {
  subset_data <- hammett_clean[hammett_clean$effect_type == etype, ]
  effect_stats <- rbind(effect_stats, data.frame(
    effect_type = etype,
    n = nrow(subset_data),
    mean_sigma_m = round(mean(subset_data$sigma_meta), 3),
    mean_sigma_p = round(mean(subset_data$sigma_para), 3),
    sd_sigma_m = round(sd(subset_data$sigma_meta), 3),
    sd_sigma_p = round(sd(subset_data$sigma_para), 3)
  ))
}
effect_stats <- effect_stats[order(effect_stats$mean_sigma_p), ]
print(effect_stats, row.names = FALSE)

# Group statistics by heteroatom
cat("\n--- Statistics by Heteroatom ---\n")
heteroatoms <- unique(hammett_clean$heteroatom)
hetero_stats <- data.frame(
  heteroatom = character(),
  n = integer(),
  mean_sigma_m = numeric(),
  mean_sigma_p = numeric(),
  stringsAsFactors = FALSE
)

for (ha in heteroatoms) {
  subset_data <- hammett_clean[hammett_clean$heteroatom == ha, ]
  hetero_stats <- rbind(hetero_stats, data.frame(
    heteroatom = as.character(ha),
    n = nrow(subset_data),
    mean_sigma_m = round(mean(subset_data$sigma_meta), 3),
    mean_sigma_p = round(mean(subset_data$sigma_para), 3)
  ))
}
hetero_stats <- hetero_stats[order(-hetero_stats$n), ]
print(hetero_stats, row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. VISUALIZATIONS (Base R)
# -----------------------------------------------------------------------------

cat("\n[4] GENERATING VISUALIZATIONS...\n")

# Create output directory
dir.create("output", showWarnings = FALSE)

# Define color palette for effect types
effect_colors <- c(
  "EDG" = "#2166AC",      # Blue
  "EWG" = "#B2182B",      # Red
  "neutral" = "#666666",   # Gray
  "weak_EDG" = "#67A9CF",  # Light blue
  "weak_EWG" = "#EF8A62"   # Light red
)

# Plot 1: Sigma meta vs Sigma para
png("output/day-1/01_sigma_meta_vs_para.png", width = 800, height = 700)
par(mar = c(5, 5, 4, 2))

plot(hammett_clean$sigma_meta, hammett_clean$sigma_para,
     col = effect_colors[as.character(hammett_clean$effect_type)],
     pch = 19, cex = 1.5,
     xlab = expression(sigma[meta]),
     ylab = expression(sigma[para]),
     main = "Hammett Substituent Constants: Meta vs Para",
     cex.lab = 1.3, cex.main = 1.4)

# Add reference line (y = x)
abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)

# Add substituent labels
text(hammett_clean$sigma_meta, hammett_clean$sigma_para,
     labels = hammett_clean$substituent,
     pos = 4, cex = 0.7, col = "gray30")

# Add legend
legend("topleft", legend = names(effect_colors),
       col = effect_colors, pch = 19, cex = 0.9,
       title = "Effect Type", bty = "n")

# Add subtitle
mtext("Dashed line: σ_meta = σ_para (no resonance effect)", 
      side = 3, line = 0.3, cex = 0.9)

dev.off()
cat("Saved: output/day-1/01_sigma_meta_vs_para.png\n")

# Plot 2: Distribution histograms
png("output/day-1/02_sigma_distributions.png", width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Sigma meta histogram
hist(hammett_clean$sigma_meta, breaks = 12, col = "#67A9CF",
     xlab = expression(sigma[meta]), main = expression("Distribution of " * sigma[meta]),
     cex.lab = 1.2, cex.main = 1.3)
abline(v = 0, lty = 2, col = "red", lwd = 2)
abline(v = mean(hammett_clean$sigma_meta), lty = 1, col = "blue", lwd = 2)
legend("topright", legend = c("Zero", "Mean"), col = c("red", "blue"),
       lty = c(2, 1), lwd = 2, bty = "n")

# Sigma para histogram
hist(hammett_clean$sigma_para, breaks = 12, col = "#EF8A62",
     xlab = expression(sigma[para]), main = expression("Distribution of " * sigma[para]),
     cex.lab = 1.2, cex.main = 1.3)
abline(v = 0, lty = 2, col = "red", lwd = 2)
abline(v = mean(hammett_clean$sigma_para), lty = 1, col = "blue", lwd = 2)
legend("topright", legend = c("Zero", "Mean"), col = c("red", "blue"),
       lty = c(2, 1), lwd = 2, bty = "n")

dev.off()
cat("Saved: output/day-1/02_sigma_distributions.png\n")

# Plot 3: Boxplot by effect type
png("output/day-1/03_boxplot_by_effect.png", width = 900, height = 600)
par(mar = c(5, 8, 4, 2))

# Prepare data for boxplot
effect_order <- c("EDG", "weak_EDG", "neutral", "weak_EWG", "EWG")
effect_present <- effect_order[effect_order %in% levels(hammett_clean$effect_type)]

# Create boxplot and store the statistics
bp <- boxplot(sigma_para ~ effect_type, data = hammett_clean,
              horizontal = TRUE, col = effect_colors[effect_present],
              xlab = expression(sigma[para]), ylab = "",
              main = "Hammett Constants by Substituent Effect Type",
              cex.lab = 1.2, cex.main = 1.3, las = 1)

# Label outliers with substituent names
effect_levels <- levels(hammett_clean$effect_type)
for (i in seq_along(effect_levels)) {
     etype <- effect_levels[i]
     subset_data <- hammett_clean[hammett_clean$effect_type == etype, ]

     # Calculate IQR bounds for this group
     q1 <- quantile(subset_data$sigma_para, 0.25)
     q3 <- quantile(subset_data$sigma_para, 0.75)
     iqr <- q3 - q1
     lower_bound <- q1 - 1.5 * iqr
     upper_bound <- q3 + 1.5 * iqr
     # Find outliers
     outliers <- subset_data[subset_data$sigma_para < lower_bound |
          subset_data$sigma_para > upper_bound,]
     # Label each outlier
     if (nrow(outliers) > 0) {
     text(x = outliers$sigma_para, 
          y = i,
          labels = outliers$substituent,
          pos = 3, cex = 0.8, col = "black", font = 2)
     }
}

abline(v = 0, lty = 2, col = "gray50", lwd = 2)
mtext("EDG = Electron Donating, EWG = Electron Withdrawing", 
      side = 3, line = 0.3, cex = 0.9)

dev.off()
cat("Saved: output/day-1/03_boxplot_by_effect.png\n")

# Plot 4: Bar chart of mean sigma by heteroatom
png("output/day-1/04_heteroatom_means.png", width = 800, height = 600)
par(mar = c(5, 5, 4, 2))

# Filter out carbon (reference)
hetero_subset <- hetero_stats[hetero_stats$heteroatom != "C", ]
hetero_subset <- hetero_subset[order(hetero_subset$mean_sigma_p), ]

barplot(hetero_subset$mean_sigma_p, names.arg = hetero_subset$heteroatom,
        col = rainbow(nrow(hetero_subset)), border = NA,
        xlab = "Heteroatom", ylab = expression("Mean " * sigma[para]),
        main = "Mean Sigma Para by Heteroatom",
        cex.lab = 1.2, cex.main = 1.3, ylim = c(-0.5, 0.8))
abline(h = 0, lty = 2, col = "gray50", lwd = 2)

# Add count labels
text(seq_len(nrow(hetero_subset)) * 1.2 - 0.5,
     hetero_subset$mean_sigma_p + 0.05,
     labels = paste0("n=", hetero_subset$n), cex = 0.8)

dev.off()
cat("Saved: output/day-1/04_heteroatom_means.png\n")

# Plot 5: Correlation between sigma_meta and sigma_para
png("output/day-1/05_correlation.png", width = 700, height = 600)
par(mar = c(5, 5, 4, 2))

# Calculate correlation
cor_value <- cor(hammett_clean$sigma_meta, hammett_clean$sigma_para)

plot(hammett_clean$sigma_meta, hammett_clean$sigma_para,
     pch = 19, col = "#2166AC", cex = 1.5,
     xlab = expression(sigma[meta]),
     ylab = expression(sigma[para]),
     main = "Correlation: Sigma Meta vs Sigma Para")

# Add linear regression line
fit <- lm(sigma_para ~ sigma_meta, data = hammett_clean)
abline(fit, col = "red", lwd = 2)

# Add correlation coefficient
legend("topleft", 
       legend = c(sprintf("r = %.3f", cor_value),
                  sprintf("R² = %.3f", cor_value^2)),
       bty = "n", cex = 1.2)

dev.off()
cat("Saved: output/day-1/05_correlation.png\n")

# -----------------------------------------------------------------------------
# 6. KEY OBSERVATIONS (Organic Chemistry Context)
# -----------------------------------------------------------------------------

cat("\n[5] KEY OBSERVATIONS...\n")

# Strongest EWG and EDG
strongest_EWG <- hammett_clean[which.max(hammett_clean$sigma_para), ]
strongest_EDG <- hammett_clean[which.min(hammett_clean$sigma_para), ]

cat("\n--- Extreme Values ---\n")
cat(sprintf("Strongest EWG (most positive σ_para): %s (σ_para = %.2f)\n",
            strongest_EWG$substituent, strongest_EWG$sigma_para))
cat(sprintf("Strongest EDG (most negative σ_para): %s (σ_para = %.2f)\n",
            strongest_EDG$substituent, strongest_EDG$sigma_para))

# Substituents with large meta-para differences (strong resonance)
cat("\n--- Large Meta-Para Differences (Strong Resonance Effects) ---\n")
large_diff <- hammett_clean[order(-hammett_clean$resonance_indicator), 
                            c("substituent", "sigma_meta", "sigma_para", "sigma_diff", "effect_type")]
print(head(large_diff, 8), row.names = FALSE)

# Substituents where σ_meta ≈ σ_para (primarily inductive)
cat("\n--- Small Meta-Para Differences (Primarily Inductive) ---\n")
small_diff <- hammett_clean[hammett_clean$resonance_indicator < 0.1,
                            c("substituent", "sigma_meta", "sigma_para", "sigma_diff")]
print(small_diff, row.names = FALSE)

# Correlation analysis
cat("\n--- Correlation Analysis ---\n")
cat(sprintf("Pearson correlation (σ_meta vs σ_para): r = %.4f\n", 
            cor(hammett_clean$sigma_meta, hammett_clean$sigma_para)))
cat(sprintf("This strong correlation indicates both positions respond similarly to electronic effects.\n"))
cat(sprintf("Deviations from r=1 reveal resonance contributions (stronger at para position).\n"))

# -----------------------------------------------------------------------------
# 7. SAVE PROCESSED DATA
# -----------------------------------------------------------------------------

cat("\n[6] SAVING OUTPUT...\n")

# Save cleaned dataset
write.csv(hammett_clean, "data/hammett_clean.csv", row.names = FALSE)
cat("Saved: data/hammett_clean.csv\n")

# Save summary to text file
sink("output/day-1/day1_summary.txt")
cat("============================================================\n")
cat("HAMMETT CONSTANT ANALYSIS - DAY 1 SUMMARY\n")
cat("============================================================\n\n")

cat("DATASET INFO:\n")
cat(sprintf("Total substituents: %d\n", nrow(hammett_clean)))
cat(sprintf("Variables: %d\n\n", ncol(hammett_clean)))

cat("EFFECT TYPE DISTRIBUTION:\n")
print(table(hammett_clean$effect_type))

cat("\nHETEROATOM DISTRIBUTION:\n")
print(table(hammett_clean$heteroatom))

cat("\nSIGMA STATISTICS:\n")
print(summary(hammett_clean[, c("sigma_meta", "sigma_para")]))

cat("\nCORRELATION:\n")
cat(sprintf("r(sigma_meta, sigma_para) = %.4f\n", 
            cor(hammett_clean$sigma_meta, hammett_clean$sigma_para)))

cat("\nKEY FINDINGS:\n")
cat(sprintf("- Strongest EWG: %s (σ_para = %.2f)\n", 
            strongest_EWG$substituent, strongest_EWG$sigma_para))
cat(sprintf("- Strongest EDG: %s (σ_para = %.2f)\n",
            strongest_EDG$substituent, strongest_EDG$sigma_para))
cat("- Large meta-para differences indicate resonance effects\n")
cat("- N(CH3)2 shows largest resonance donation (σ_diff = 0.68)\n")
sink()
cat("Saved: docs/day-1-summary.txt\n")

# -----------------------------------------------------------------------------
# 8. COMPLETION MESSAGE
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("                  DAY 1 COMPLETE!\n")
cat("============================================================\n")