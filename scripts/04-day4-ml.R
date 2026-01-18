# =============================================================================
# Hammett Substituent Constant Analysis
# Day 4: Machine Learning & Predictive Modeling
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("   HAMMETT SUBSTITUENT CONSTANT ANALYSIS - DAY 4\n")
cat("   Machine Learning & Predictive Modeling\n")
cat("============================================================\n")

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("\n[1] LOADING DATA...\n")

hammett <- read.csv("data/hammett_day_3_r.csv", stringsAsFactors = FALSE)

# Also load extended dataset for more features
hammett_ext <- read.csv("data/hammett_extended.csv", stringsAsFactors = FALSE)

cat(sprintf("Loaded %d substituents from Day 3 dataset\n", nrow(hammett)))

# Merge datasets to get sigma_plus and sigma_minus
hammett_merged <- merge(hammett,
                        hammett_ext[, c("substituent", "sigma_plus", "sigma_minus")],
                        by = "substituent", all.x = TRUE)

cat(sprintf("After merge: %d substituents with extended features\n", nrow(hammett_merged)))

# -----------------------------------------------------------------------------
# 2. FEATURE ENGINEERING
# -----------------------------------------------------------------------------

cat("\n[2] FEATURE ENGINEERING...\n")

# Create additional derived features
hammett_merged$sigma_meta_sq <- hammett_merged$sigma_meta^2  # Quadratic term
hammett_merged$sigma_I_sq <- hammett_merged$sigma_I^2
hammett_merged$meta_I_interaction <- hammett_merged$sigma_meta * hammett_merged$sigma_I

# Binary features
hammett_merged$is_EDG <- as.integer(hammett_merged$effect_type == "EDG")
hammett_merged$is_EWG <- as.integer(hammett_merged$effect_type == "EWG")
hammett_merged$has_N <- as.integer(hammett_merged$heteroatom == "N")
hammett_merged$has_O <- as.integer(hammett_merged$heteroatom == "O")
hammett_merged$is_halogen <- as.integer(hammett_merged$heteroatom %in% c("F", "Cl", "Br", "I"))

cat("Created features:\n")
cat("  - sigma_meta_sq (quadratic term)\n")
cat("  - sigma_I_sq (quadratic term)\n")
cat("  - meta_I_interaction (interaction term)\n")
cat("  - Binary indicators: is_EDG, is_EWG, has_N, has_O, is_halogen\n")

# -----------------------------------------------------------------------------
# 3. TRAIN/TEST SPLIT
# -----------------------------------------------------------------------------

cat("\n[3] TRAIN/TEST SPLIT...\n")

set.seed(42)  # For reproducibility

# 80/20 split
n <- nrow(hammett_merged)
train_idx <- sample(1:n, size = floor(0.8 * n))
test_idx <- setdiff(1:n, train_idx)

train_data <- hammett_merged[train_idx, ]
test_data <- hammett_merged[test_idx, ]

cat(sprintf("Training set: %d samples (%.0f%%)\n", nrow(train_data), 100 * nrow(train_data) / n))
cat(sprintf("Test set: %d samples (%.0f%%)\n", nrow(test_data), 100 * nrow(test_data) / n))

# Visualize the split
png("output/day-4/15_train_test_split.png", width = 900, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Plot 1: Distribution comparison
hist(train_data$sigma_para, breaks = 10, col = "#2166AC80", border = "white",
     main = "Train vs Test Distribution",
     xlab = expression(sigma[para]), xlim = range(hammett_merged$sigma_para),
     ylim = c(0, 12))
hist(test_data$sigma_para, breaks = 8, col = "#B2182B80", border = "white", add = TRUE)
legend("topright", legend = c("Train", "Test"), fill = c("#2166AC80", "#B2182B80"), bty = "n")

# Plot 2: Scatter showing split
plot(hammett_merged$sigma_meta, hammett_merged$sigma_para,
     pch = 19, cex = 1.5,
     col = ifelse(1:n %in% train_idx, "#2166AC", "#B2182B"),
     xlab = expression(sigma[meta]),
     ylab = expression(sigma[para]),
     main = "Train/Test Split Visualization")
legend("topleft", legend = c("Train", "Test"), col = c("#2166AC", "#B2182B"),
       pch = 19, bty = "n")

dev.off()
cat("Saved: output/day-4/15_train_test_split.png\n")

# -----------------------------------------------------------------------------
# 4. MODEL 1: SIMPLE LINEAR REGRESSION (Baseline)
# -----------------------------------------------------------------------------

cat("\n[4] MODEL 1: SIMPLE LINEAR REGRESSION (Baseline)...\n")
cat("─────────────────────────────────────────────────────────────\n")

model1 <- lm(sigma_para ~ sigma_meta, data = train_data)

cat("\n--- Model 1: σ_para ~ σ_meta ---\n")
print(summary(model1))

# Predictions
train_data$pred_m1 <- predict(model1, train_data)
test_data$pred_m1 <- predict(model1, test_data)

# Calculate metrics
calc_metrics <- function(actual, predicted) {
  residuals <- actual - predicted
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(residuals))
  ss_res <- sum(residuals^2)
  ss_tot <- sum((actual - mean(actual))^2)
  r_squared <- 1 - ss_res / ss_tot
  
  return(list(MSE = mse, RMSE = rmse, MAE = mae, R2 = r_squared))
}

metrics_m1_train <- calc_metrics(train_data$sigma_para, train_data$pred_m1)
metrics_m1_test <- calc_metrics(test_data$sigma_para, test_data$pred_m1)

cat("\n--- Model 1 Performance ---\n")
cat(sprintf("Training: R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m1_train$R2, metrics_m1_train$RMSE, metrics_m1_train$MAE))
cat(sprintf("Test:     R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m1_test$R2, metrics_m1_test$RMSE, metrics_m1_test$MAE))

# Visualization
png("output/day-4/16_model1_baseline.png", width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Training fit
plot(train_data$sigma_para, train_data$pred_m1,
     pch = 19, col = "#2166AC", cex = 1.5,
     xlab = expression("Actual " * sigma[para]),
     ylab = expression("Predicted " * sigma[para]),
     main = "Model 1: Training Set",
     xlim = range(hammett_merged$sigma_para),
     ylim = range(hammett_merged$sigma_para))
abline(0, 1, col = "#B2182B", lwd = 2, lty = 2)
text(-0.6, 0.7, sprintf("R² = %.3f\nRMSE = %.3f", 
                         metrics_m1_train$R2, metrics_m1_train$RMSE),
     cex = 1.1, adj = 0)

# Test fit
plot(test_data$sigma_para, test_data$pred_m1,
     pch = 19, col = "#B2182B", cex = 1.5,
     xlab = expression("Actual " * sigma[para]),
     ylab = expression("Predicted " * sigma[para]),
     main = "Model 1: Test Set",
     xlim = range(hammett_merged$sigma_para),
     ylim = range(hammett_merged$sigma_para))
abline(0, 1, col = "#2166AC", lwd = 2, lty = 2)
text(-0.6, 0.7, sprintf("R² = %.3f\nRMSE = %.3f",
                         metrics_m1_test$R2, metrics_m1_test$RMSE),
     cex = 1.1, adj = 0)

# Label test points
text(test_data$sigma_para, test_data$pred_m1, 
     labels = test_data$substituent, pos = 4, cex = 0.7)

dev.off()
cat("Saved: output/day-4/16_model1_baseline.png\n")

# -----------------------------------------------------------------------------
# 5. MODEL 2: MULTIPLE LINEAR REGRESSION
# -----------------------------------------------------------------------------

cat("\n[5] MODEL 2: MULTIPLE LINEAR REGRESSION...\n")
cat("─────────────────────────────────────────────────────────────\n")

# Use sigma_meta and sigma_I
model2 <- lm(sigma_para ~ sigma_meta + sigma_I, data = train_data)

cat("\n--- Model 2: σ_para ~ σ_meta + σ_I ---\n")
print(summary(model2))

# Predictions
train_data$pred_m2 <- predict(model2, train_data)
test_data$pred_m2 <- predict(model2, test_data)

metrics_m2_train <- calc_metrics(train_data$sigma_para, train_data$pred_m2)
metrics_m2_test <- calc_metrics(test_data$sigma_para, test_data$pred_m2)

cat("\n--- Model 2 Performance ---\n")
cat(sprintf("Training: R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m2_train$R2, metrics_m2_train$RMSE, metrics_m2_train$MAE))
cat(sprintf("Test:     R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m2_test$R2, metrics_m2_test$RMSE, metrics_m2_test$MAE))

# Visualization: Coefficient plot
png("output/day-4/17_model2_coefficients.png", width = 800, height = 500)
par(mar = c(5, 8, 4, 2))

coef_m2 <- summary(model2)$coefficients
coef_names <- c("Intercept", expression(sigma[meta]), expression(sigma[I]))
coef_vals <- coef_m2[, 1]
coef_se <- coef_m2[, 2]

# Horizontal bar plot
bp <- barplot(coef_vals, names.arg = c("Intercept", "sigma_meta", "sigma_resonance"),
              horiz = TRUE, col = c("#999999", "#2166AC", "#B2182B"),
              main = "Model 2: Coefficient Estimates with 95% CI",
              xlab = "Coefficient Value", xlim = c(-2, 3),
              las = 1, border = NA)

# Add error bars (95% CI = ±1.96*SE)
arrows(coef_vals - 1.96 * coef_se, bp, coef_vals + 1.96 * coef_se, bp,
       angle = 90, code = 3, length = 0.1, lwd = 2)

# Add zero reference line
abline(v = 0, lty = 2, col = "gray50")

# Add coefficient values as text
text(coef_vals + 0.3, bp, sprintf("%.3f", coef_vals), cex = 0.9)

dev.off()
cat("Saved: output/day-4/17_model2_coefficients.png\n")

# Actual vs Predicted for Model 2
png("output/day-4/18_model2_predictions.png", width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

plot(train_data$sigma_para, train_data$pred_m2,
     pch = 19, col = "#2166AC", cex = 1.5,
     xlab = expression("Actual " * sigma[para]),
     ylab = expression("Predicted " * sigma[para]),
     main = "Model 2: Training Set",
     xlim = range(hammett_merged$sigma_para),
     ylim = range(hammett_merged$sigma_para))
abline(0, 1, col = "#B2182B", lwd = 2, lty = 2)
text(-0.6, 0.7, sprintf("R² = %.3f\nRMSE = %.3f",
                         metrics_m2_train$R2, metrics_m2_train$RMSE),
     cex = 1.1, adj = 0)

plot(test_data$sigma_para, test_data$pred_m2,
     pch = 19, col = "#B2182B", cex = 1.5,
     xlab = expression("Actual " * sigma[para]),
     ylab = expression("Predicted " * sigma[para]),
     main = "Model 2: Test Set",
     xlim = range(hammett_merged$sigma_para),
     ylim = range(hammett_merged$sigma_para))
abline(0, 1, col = "#2166AC", lwd = 2, lty = 2)
text(-0.6, 0.7, sprintf("R² = %.3f\nRMSE = %.3f",
                         metrics_m2_test$R2, metrics_m2_test$RMSE),
     cex = 1.1, adj = 0)
text(test_data$sigma_para, test_data$pred_m2,
     labels = test_data$substituent, pos = 4, cex = 0.7)

dev.off()
cat("Saved: output/day-4/18_model2_predictions.png\n")

# -----------------------------------------------------------------------------
# 6. MODEL 3: POLYNOMIAL REGRESSION (Quadratic)
# -----------------------------------------------------------------------------

cat("\n[6] MODEL 3: POLYNOMIAL REGRESSION...\n")
cat("─────────────────────────────────────────────────────────────\n")

model3 <- lm(sigma_para ~ sigma_meta + sigma_meta_sq, data = train_data)

cat("\n--- Model 3: σ_para ~ σ_meta + σ_meta² ---\n")
print(summary(model3))

train_data$pred_m3 <- predict(model3, train_data)
test_data$pred_m3 <- predict(model3, test_data)

metrics_m3_train <- calc_metrics(train_data$sigma_para, train_data$pred_m3)
metrics_m3_test <- calc_metrics(test_data$sigma_para, test_data$pred_m3)

cat("\n--- Model 3 Performance ---\n")
cat(sprintf("Training: R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m3_train$R2, metrics_m3_train$RMSE, metrics_m3_train$MAE))
cat(sprintf("Test:     R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m3_test$R2, metrics_m3_test$RMSE, metrics_m3_test$MAE))

# Visualization: Polynomial fit curve
png("output/day-4/19_model3_polynomial.png", width = 900, height = 600)
par(mar = c(5, 5, 4, 2))

plot(hammett_merged$sigma_meta, hammett_merged$sigma_para,
     pch = 19, cex = 1.5,
     col = ifelse(1:n %in% train_idx, "#2166AC", "#B2182B"),
     xlab = expression(sigma[meta]),
     ylab = expression(sigma[para]),
     main = "Model 3: Polynomial Regression (Quadratic)")

# Add linear fit (Model 1)
abline(model1, col = "#67A9CF", lwd = 2, lty = 2)

# Add polynomial fit curve
meta_seq <- seq(min(hammett_merged$sigma_meta), max(hammett_merged$sigma_meta), length.out = 100)
poly_pred <- predict(model3, newdata = data.frame(sigma_meta = meta_seq, 
                                                   sigma_meta_sq = meta_seq^2))
lines(meta_seq, poly_pred, col = "#B2182B", lwd = 3)

# Add reference line y = x
abline(0, 1, col = "gray50", lty = 3)

legend("topleft", 
       legend = c("Train data", "Test data", "Linear fit", "Quadratic fit", "y = x"),
       col = c("#2166AC", "#B2182B", "#67A9CF", "#B2182B", "gray50"),
       pch = c(19, 19, NA, NA, NA), lty = c(NA, NA, 2, 1, 3), lwd = 2,
       bty = "n")

# Add equation
coef3 <- coef(model3)
eq_text <- sprintf("y = %.3f + %.3f·x + %.3f·x²", coef3[1], coef3[2], coef3[3])
text(0.4, -0.5, eq_text, cex = 1.1)

dev.off()
cat("Saved: output/day-4/19_model3_polynomial.png\n")

# -----------------------------------------------------------------------------
# 7. MODEL 4: MULTIPLE REGRESSION WITH INTERACTION
# -----------------------------------------------------------------------------

cat("\n[7] MODEL 4: REGRESSION WITH INTERACTION TERM...\n")
cat("─────────────────────────────────────────────────────────────\n")

model4 <- lm(sigma_para ~ sigma_meta * sigma_I, data = train_data)

cat("\n--- Model 4: σ_para ~ σ_meta * σ_I (with interaction) ---\n")
print(summary(model4))

train_data$pred_m4 <- predict(model4, train_data)
test_data$pred_m4 <- predict(model4, test_data)

metrics_m4_train <- calc_metrics(train_data$sigma_para, train_data$pred_m4)
metrics_m4_test <- calc_metrics(test_data$sigma_para, test_data$pred_m4)

cat("\n--- Model 4 Performance ---\n")
cat(sprintf("Training: R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m4_train$R2, metrics_m4_train$RMSE, metrics_m4_train$MAE))
cat(sprintf("Test:     R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m4_test$R2, metrics_m4_test$RMSE, metrics_m4_test$MAE))

# Visualization: Interaction effect
png("output/day-4/20_model4_interaction.png", width = 900, height = 600)
par(mar = c(5, 5, 4, 2))

# Create grid for 3D-like visualization
# Show predicted sigma_para at different sigma_I levels

sigma_I_levels <- quantile(train_data$sigma_I, c(0.25, 0.5, 0.75), na.rm = TRUE)
colors <- c("#2166AC", "#666666", "#B2182B")

plot(NULL, xlim = range(hammett_merged$sigma_meta, na.rm = TRUE),
     ylim = range(hammett_merged$sigma_para, na.rm = TRUE),
     xlab = expression(sigma[meta]),
     ylab = expression("Predicted " * sigma[para]),
     main = "Model 4: Interaction Effect (σ_meta × σ_I)")

meta_seq <- seq(min(hammett_merged$sigma_meta), max(hammett_merged$sigma_meta), length.out = 100)

for (i in 1:3) {
  pred_vals <- predict(model4, newdata = data.frame(
    sigma_meta = meta_seq,
    sigma_I = rep(sigma_I_levels[i], 100)
  ))
  lines(meta_seq, pred_vals, col = colors[i], lwd = 3)
}

# Add actual data points
points(hammett_merged$sigma_meta, hammett_merged$sigma_para,
       pch = 19, col = "#00000050", cex = 1.2)

legend("topleft",
       legend = c(sprintf("σ_I = %.2f (25th %%ile)", sigma_I_levels[1]),
                  sprintf("σ_I = %.2f (median)", sigma_I_levels[2]),
                  sprintf("σ_I = %.2f (75th %%ile)", sigma_I_levels[3]),
                  "Actual data"),
       col = c(colors, "#00000050"), lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19),
       lwd = 2, bty = "n")

dev.off()
cat("Saved: output/day-4/20_model4_interaction.png\n")

# -----------------------------------------------------------------------------
# 8. MODEL 5: KITCHEN SINK MODEL (All features)
# -----------------------------------------------------------------------------

cat("\n[8] MODEL 5: FULL MODEL (Multiple predictors)...\n")
cat("─────────────────────────────────────────────────────────────\n")

# Use all available numeric predictors
model5 <- lm(sigma_para ~ sigma_meta + sigma_I + pi + MR + is_EDG + is_EWG, 
             data = train_data)

cat("\n--- Model 5: Full model with multiple predictors ---\n")
print(summary(model5))

train_data$pred_m5 <- predict(model5, train_data)
test_data$pred_m5 <- predict(model5, test_data)

metrics_m5_train <- calc_metrics(train_data$sigma_para, train_data$pred_m5)
metrics_m5_test <- calc_metrics(test_data$sigma_para, test_data$pred_m5)

cat("\n--- Model 5 Performance ---\n")
cat(sprintf("Training: R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m5_train$R2, metrics_m5_train$RMSE, metrics_m5_train$MAE))
cat(sprintf("Test:     R² = %.4f, RMSE = %.4f, MAE = %.4f\n",
            metrics_m5_test$R2, metrics_m5_test$RMSE, metrics_m5_test$MAE))

# Variable importance (standardized coefficients)
png("output/day-4/21_model5_importance.png", width = 900, height = 600)
par(mar = c(5, 10, 4, 2))

# Calculate standardized coefficients (beta weights)
# Beta = b * (SD_x / SD_y)
coef5 <- coef(model5)[-1]  # Exclude intercept
predictors <- c("sigma_meta", "sigma_I", "pi", "MR", "is_EDG", "is_EWG")
sd_y <- sd(train_data$sigma_para)

std_coefs <- sapply(seq_along(predictors), function(i) {
  coef5[i] * sd(train_data[[predictors[i]]], na.rm = TRUE) / sd_y
})
names(std_coefs) <- predictors

# Sort by absolute value
std_coefs_sorted <- std_coefs[order(abs(std_coefs))]

barplot(std_coefs_sorted, horiz = TRUE, las = 1,
        col = ifelse(std_coefs_sorted > 0, "#B2182B", "#2166AC"),
        main = "Model 5: Standardized Coefficients (Variable Importance)",
        xlab = "Standardized Coefficient (β)",
        border = NA)
abline(v = 0, lty = 2)

dev.off()
cat("Saved: output/day-4/21_model5_importance.png\n")

# Model 5 Predictions Plot
set.seed(42)
n <- nrow(hammett)
train_idx <- sample(1:n, size = floor(0.8 * n))
test_idx <- setdiff(1:n, train_idx)

train_data <- hammett[train_idx, ]
test_data <- hammett[test_idx, ]

# Create binary features
train_data$is_EDG <- as.integer(train_data$effect_type == "EDG")
train_data$is_EWG <- as.integer(train_data$effect_type == "EWG")
test_data$is_EDG <- as.integer(test_data$effect_type == "EDG")
test_data$is_EWG <- as.integer(test_data$effect_type == "EWG")

# Fit Model 5
model5 <- lm(sigma_para ~ sigma_meta + sigma_I + pi + MR + is_EDG + is_EWG, 
             data = train_data)

# Predictions
train_data$pred_m5 <- predict(model5, train_data)
test_data$pred_m5 <- predict(model5, test_data)

# Calculate metrics
calc_metrics <- function(actual, predicted) {
  residuals <- actual - predicted
  ss_res <- sum(residuals^2)
  ss_tot <- sum((actual - mean(actual))^2)
  r_squared <- 1 - ss_res / ss_tot
  rmse <- sqrt(mean(residuals^2))
  return(list(R2 = r_squared, RMSE = rmse))
}

metrics_train <- calc_metrics(train_data$sigma_para, train_data$pred_m5)
metrics_test <- calc_metrics(test_data$sigma_para, test_data$pred_m5)

# Generate plot
png("output/day-4/22_model5_predictions.png", width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Training set
plot(train_data$sigma_para, train_data$pred_m5,
     pch = 19, col = "#2166AC", cex = 1.5,
     xlab = expression("Actual " * sigma[para]),
     ylab = expression("Predicted " * sigma[para]),
     main = "Model 5: Training Set",
     xlim = range(hammett$sigma_para),
     ylim = range(hammett$sigma_para))
abline(0, 1, col = "#B2182B", lwd = 2, lty = 2)
text(-0.6, 0.7, sprintf("R² = %.3f\nRMSE = %.3f", 
                         metrics_train$R2, metrics_train$RMSE),
     cex = 1.1, adj = 0)

# Test set
plot(test_data$sigma_para, test_data$pred_m5,
     pch = 19, col = "#B2182B", cex = 1.5,
     xlab = expression("Actual " * sigma[para]),
     ylab = expression("Predicted " * sigma[para]),
     main = "Model 5: Test Set",
     xlim = range(hammett$sigma_para),
     ylim = range(hammett$sigma_para))
abline(0, 1, col = "#2166AC", lwd = 2, lty = 2)
text(-0.6, 0.7, sprintf("R² = %.3f\nRMSE = %.3f",
                         metrics_test$R2, metrics_test$RMSE),
     cex = 1.1, adj = 0)
text(test_data$sigma_para, test_data$pred_m5,
     labels = test_data$substituent, pos = 4, cex = 0.7)

dev.off()
cat("Saved: output/day-4/22_model5_predictions.png\n")

# -----------------------------------------------------------------------------
# 9. MODEL COMPARISON
# -----------------------------------------------------------------------------

cat("\n[9] MODEL COMPARISON...\n")
cat("─────────────────────────────────────────────────────────────\n")

# Compile all results
model_comparison <- data.frame(
  Model = c("M1: Simple (σ_meta)",
            "M2: Multiple (σ_meta + σ_I)",
            "M3: Polynomial (σ_meta²)",
            "M4: Interaction (σ_meta × σ_I)",
            "M5: Full model"),
  Train_R2 = c(metrics_m1_train$R2, metrics_m2_train$R2, metrics_m3_train$R2,
               metrics_m4_train$R2, metrics_m5_train$R2),
  Test_R2 = c(metrics_m1_test$R2, metrics_m2_test$R2, metrics_m3_test$R2,
              metrics_m4_test$R2, metrics_m5_test$R2),
  Train_RMSE = c(metrics_m1_train$RMSE, metrics_m2_train$RMSE, metrics_m3_train$RMSE,
                 metrics_m4_train$RMSE, metrics_m5_train$RMSE),
  Test_RMSE = c(metrics_m1_test$RMSE, metrics_m2_test$RMSE, metrics_m3_test$RMSE,
                metrics_m4_test$RMSE, metrics_m5_test$RMSE),
  Num_Predictors = c(1, 2, 2, 3, 6)
)

cat("\n--- Model Comparison Summary ---\n")
print(model_comparison, row.names = FALSE)

# Calculate overfitting measure (Train R² - Test R²)
model_comparison$Overfit <- model_comparison$Train_R2 - model_comparison$Test_R2

# Visualization: Model comparison
png("output/day-4/23_model_comparison.png", width = 1000, height = 600)
par(mfrow = c(1, 2), mar = c(8, 5, 4, 2))

# Plot 1: R² comparison
bar_positions <- barplot(rbind(model_comparison$Train_R2, model_comparison$Test_R2),
                         beside = TRUE, col = c("#2166AC", "#B2182B"),
                         names.arg = c("M1", "M2", "M3", "M4", "M5"),
                         ylim = c(0, 1.1),
                         main = "Model Comparison: R²",
                         ylab = "R²")
legend("topright", legend = c("Train", "Test"), fill = c("#2166AC", "#B2182B"), bty = "n")

# Add value labels
text(bar_positions[1, ], model_comparison$Train_R2 + 0.05, 
     sprintf("%.2f", model_comparison$Train_R2), cex = 0.8)
text(bar_positions[2, ], model_comparison$Test_R2 + 0.05,
     sprintf("%.2f", model_comparison$Test_R2), cex = 0.8)

# Plot 2: RMSE comparison
bar_positions2 <- barplot(rbind(model_comparison$Train_RMSE, model_comparison$Test_RMSE),
                          beside = TRUE, col = c("#2166AC", "#B2182B"),
                          names.arg = c("M1", "M2", "M3", "M4", "M5"),
                          ylim = c(0, max(model_comparison$Test_RMSE) * 1.3),
                          main = "Model Comparison: RMSE",
                          ylab = "RMSE")
legend("topright", legend = c("Train", "Test"), fill = c("#2166AC", "#B2182B"), bty = "n")

dev.off()
cat("Saved: output/day-4/23_model_comparison.png\n")

# -----------------------------------------------------------------------------
# 10. CROSS-VALIDATION
# -----------------------------------------------------------------------------

cat("\n[10] K-FOLD CROSS-VALIDATION...\n")
cat("─────────────────────────────────────────────────────────────\n")

# Implement k-fold CV manually (no external packages)
k_fold_cv <- function(data, formula, k = 5) {
  set.seed(123)
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  
  cv_results <- data.frame(
    Fold = 1:k,
    Train_R2 = numeric(k),
    Test_R2 = numeric(k),
    Test_RMSE = numeric(k)
  )
  
  for (i in 1:k) {
    train_fold <- data[folds != i, ]
    test_fold <- data[folds == i, ]
    
    model <- lm(formula, data = train_fold)
    
    train_pred <- predict(model, train_fold)
    test_pred <- predict(model, test_fold)
    
    cv_results$Train_R2[i] <- 1 - sum((train_fold$sigma_para - train_pred)^2) / 
                                  sum((train_fold$sigma_para - mean(train_fold$sigma_para))^2)
    cv_results$Test_R2[i] <- 1 - sum((test_fold$sigma_para - test_pred)^2) /
                                 sum((test_fold$sigma_para - mean(test_fold$sigma_para))^2)
    cv_results$Test_RMSE[i] <- sqrt(mean((test_fold$sigma_para - test_pred)^2))
  }
  
  return(cv_results)
}

cat("\n--- 5-Fold Cross-Validation: Model 2 (σ_meta + σ_I) ---\n")
cv_results_m2 <- k_fold_cv(hammett_merged, sigma_para ~ sigma_meta + sigma_I, k = 5)
print(cv_results_m2)
cat(sprintf("\nMean Test R²: %.4f (SD = %.4f)\n", 
            mean(cv_results_m2$Test_R2), sd(cv_results_m2$Test_R2)))
cat(sprintf("Mean Test RMSE: %.4f (SD = %.4f)\n",
            mean(cv_results_m2$Test_RMSE), sd(cv_results_m2$Test_RMSE)))

# CV for all models
cv_m1 <- k_fold_cv(hammett_merged, sigma_para ~ sigma_meta, k = 5)
cv_m2 <- k_fold_cv(hammett_merged, sigma_para ~ sigma_meta + sigma_I, k = 5)
cv_m3 <- k_fold_cv(hammett_merged, sigma_para ~ sigma_meta + I(sigma_meta^2), k = 5)

cv_summary <- data.frame(
  Model = c("M1: Simple", "M2: Multiple", "M3: Polynomial"),
  Mean_CV_R2 = c(mean(cv_m1$Test_R2), mean(cv_m2$Test_R2), mean(cv_m3$Test_R2)),
  SD_CV_R2 = c(sd(cv_m1$Test_R2), sd(cv_m2$Test_R2), sd(cv_m3$Test_R2)),
  Mean_CV_RMSE = c(mean(cv_m1$Test_RMSE), mean(cv_m2$Test_RMSE), mean(cv_m3$Test_RMSE))
)

cat("\n--- Cross-Validation Summary ---\n")
print(cv_summary, row.names = FALSE)

# Visualization: CV results
png("output/day-4/24_cross_validation.png", width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Plot 1: CV R² by fold
plot(1:5, cv_m1$Test_R2, type = "b", pch = 19, col = "#2166AC", lwd = 2,
     ylim = c(0.5, 1), xlab = "Fold", ylab = "Test R²",
     main = "5-Fold CV: Test R² by Fold")
lines(1:5, cv_m2$Test_R2, type = "b", pch = 17, col = "#B2182B", lwd = 2)
lines(1:5, cv_m3$Test_R2, type = "b", pch = 15, col = "#4DAF4A", lwd = 2)
legend("bottomright", legend = c("M1: Simple", "M2: Multiple", "M3: Polynomial"),
       col = c("#2166AC", "#B2182B", "#4DAF4A"), pch = c(19, 17, 15), lwd = 2, bty = "n")

# Plot 2: Mean CV R² with error bars
bp <- barplot(cv_summary$Mean_CV_R2, names.arg = c("M1", "M2", "M3"),
              col = c("#2166AC", "#B2182B", "#4DAF4A"),
              ylim = c(0, 1.1), main = "Mean CV R² (±1 SD)",
              ylab = "Mean Test R²")
arrows(bp, cv_summary$Mean_CV_R2 - cv_summary$SD_CV_R2,
       bp, cv_summary$Mean_CV_R2 + cv_summary$SD_CV_R2,
       angle = 90, code = 3, length = 0.1, lwd = 2)
text(bp, cv_summary$Mean_CV_R2 + cv_summary$SD_CV_R2 + 0.05,
     sprintf("%.3f", cv_summary$Mean_CV_R2), cex = 0.9)

dev.off()
cat("Saved: output/day-4/24_cross_validation.png\n")

# -----------------------------------------------------------------------------
# 11. RESIDUAL DIAGNOSTICS (Best Model)
# -----------------------------------------------------------------------------

cat("\n[11] RESIDUAL DIAGNOSTICS FOR BEST MODEL...\n")
cat("─────────────────────────────────────────────────────────────\n")

# Use Model 2 as best balance of performance and simplicity
best_model <- model2
hammett_merged$best_pred <- predict(best_model, hammett_merged)
hammett_merged$best_resid <- hammett_merged$sigma_para - hammett_merged$best_pred

png("output/day-4/25_residual_diagnostics.png", width = 1000, height = 800)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

# Plot 1: Residuals vs Fitted
plot(hammett_merged$best_pred, hammett_merged$best_resid,
     pch = 19, col = "#2166AC", cex = 1.2,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red", lwd = 2)
lines(lowess(hammett_merged$best_pred, hammett_merged$best_resid), col = "blue", lwd = 2)

# Plot 2: Q-Q Plot
qqnorm(hammett_merged$best_resid, pch = 19, col = "#2166AC",
       main = "Normal Q-Q Plot of Residuals")
qqline(hammett_merged$best_resid, col = "red", lwd = 2)

# Plot 3: Scale-Location (sqrt of standardized residuals)
std_resid <- hammett_merged$best_resid / sd(hammett_merged$best_resid)
plot(hammett_merged$best_pred, sqrt(abs(std_resid)),
     pch = 19, col = "#2166AC", cex = 1.2,
     xlab = "Fitted Values", ylab = expression(sqrt("|Standardized Residuals|")),
     main = "Scale-Location Plot")
lines(lowess(hammett_merged$best_pred, sqrt(abs(std_resid))), col = "blue", lwd = 2)

# Plot 4: Residuals by substituent
resid_order <- order(hammett_merged$best_resid)
barplot(hammett_merged$best_resid[resid_order],
        names.arg = hammett_merged$substituent[resid_order],
        col = ifelse(hammett_merged$best_resid[resid_order] > 0, "#EF8A62", "#67A9CF"),
        las = 2, cex.names = 0.5,
        main = "Residuals by Substituent",
        ylab = "Residual")
abline(h = 0)

dev.off()
cat("Saved: output/day-4/25_residual_diagnostics.png\n")

# Identify influential points
cat("\n--- Influential Points (|Residual| > 2*SD) ---\n")
threshold <- 2 * sd(hammett_merged$best_resid)
influential <- hammett_merged[abs(hammett_merged$best_resid) > threshold,
                               c("substituent", "sigma_para", "best_pred", "best_resid")]
if (nrow(influential) > 0) {
  print(influential, row.names = FALSE)
} else {
  cat("No points exceed 2*SD threshold\n")
}

# -----------------------------------------------------------------------------
# 12. PREDICTION INTERVALS
# -----------------------------------------------------------------------------

cat("\n[12] PREDICTION INTERVALS...\n")
cat("─────────────────────────────────────────────────────────────\n")

# Get prediction intervals for test set
pred_intervals <- predict(best_model, test_data, interval = "prediction", level = 0.95)
test_data$pred_lower <- pred_intervals[, "lwr"]
test_data$pred_upper <- pred_intervals[, "upr"]

cat("\n--- 95% Prediction Intervals for Test Set ---\n")
print(test_data[, c("substituent", "sigma_para", "pred_m2", "pred_lower", "pred_upper")],
      row.names = FALSE)

# Check coverage
coverage <- mean(test_data$sigma_para >= test_data$pred_lower & 
                  test_data$sigma_para <= test_data$pred_upper)
cat(sprintf("\nActual coverage: %.1f%% (expected: 95%%)\n", coverage * 100))

# Visualization
png("output/day-4/26_prediction_intervals.png", width = 900, height = 600)
par(mar = c(8, 5, 4, 2))

# Order by predicted value
test_order <- order(test_data$pred_m2)
test_ordered <- test_data[test_order, ]

x_pos <- seq_len(nrow(test_ordered))

plot(x_pos, test_ordered$pred_m2, pch = 19, col = "#2166AC", cex = 1.5,
     ylim = range(c(test_ordered$pred_lower, test_ordered$pred_upper, test_ordered$sigma_para)),
     xlab = "", ylab = expression(sigma[para]),
     main = "95% Prediction Intervals for Test Set",
     xaxt = "n")

# Add prediction intervals
arrows(x_pos, test_ordered$pred_lower, x_pos, test_ordered$pred_upper,
       angle = 90, code = 3, length = 0.1, col = "#2166AC", lwd = 2)

# Add actual values
points(x_pos, test_ordered$sigma_para, pch = 17, col = "#B2182B", cex = 1.5)

# X-axis labels
axis(1, at = x_pos, labels = test_ordered$substituent, las = 2, cex.axis = 0.8)

legend("topleft", legend = c("Predicted", "Actual", "95% PI"),
       col = c("#2166AC", "#B2182B", "#2166AC"),
       pch = c(19, 17, NA), lty = c(NA, NA, 1), lwd = 2, bty = "n")

dev.off()
cat("Saved: output/day-4/26_prediction_intervals.png\n")

# -----------------------------------------------------------------------------
# 13. LEARNING CURVE
# -----------------------------------------------------------------------------

cat("\n[13] LEARNING CURVE ANALYSIS...\n")
cat("─────────────────────────────────────────────────────────────\n")

# How does performance change with training set size?
learning_curve <- function(data, formula, train_sizes, n_iter = 10) {
  results <- data.frame(
    Train_Size = train_sizes,
    Mean_Train_R2 = numeric(length(train_sizes)),
    Mean_Test_R2 = numeric(length(train_sizes)),
    SD_Test_R2 = numeric(length(train_sizes))
  )
  
  for (i in seq_along(train_sizes)) {
    size <- train_sizes[i]
    train_r2 <- numeric(n_iter)
    test_r2 <- numeric(n_iter)
    
    for (j in 1:n_iter) {
      # Random split
      idx <- sample(seq_len(nrow(data)), size)
      train <- data[idx, ]
      test <- data[-idx, ]
      
      if (nrow(test) < 3) next
      
      model <- lm(formula, data = train)
      
      train_pred <- predict(model, train)
      test_pred <- predict(model, test)
      
      train_r2[j] <- 1 - sum((train$sigma_para - train_pred)^2) /
                         sum((train$sigma_para - mean(train$sigma_para))^2)
      test_r2[j] <- 1 - sum((test$sigma_para - test_pred)^2) /
                        sum((test$sigma_para - mean(test$sigma_para))^2)
    }
    
    results$Mean_Train_R2[i] <- mean(train_r2, na.rm = TRUE)
    results$Mean_Test_R2[i] <- mean(test_r2, na.rm = TRUE)
    results$SD_Test_R2[i] <- sd(test_r2, na.rm = TRUE)
  }
  
  return(results)
}

train_sizes <- c(8, 12, 16, 20, 24, 28)
lc_results <- learning_curve(hammett_merged, sigma_para ~ sigma_meta + sigma_I, 
                              train_sizes, n_iter = 20)

cat("\n--- Learning Curve Results ---\n")
print(lc_results, row.names = FALSE)

png("output/day-4/26_learning_curve.png", width = 800, height = 600)
par(mar = c(5, 5, 4, 2))

plot(lc_results$Train_Size, lc_results$Mean_Train_R2, type = "b",
     pch = 19, col = "#2166AC", lwd = 2, ylim = c(0.5, 1),
     xlab = "Training Set Size", ylab = "R²",
     main = "Learning Curve: Model 2 (σ_meta + σ_I)")

lines(lc_results$Train_Size, lc_results$Mean_Test_R2, type = "b",
      pch = 17, col = "#B2182B", lwd = 2)

# Add error bands for test
polygon(c(lc_results$Train_Size, rev(lc_results$Train_Size)),
        c(lc_results$Mean_Test_R2 + lc_results$SD_Test_R2,
          rev(lc_results$Mean_Test_R2 - lc_results$SD_Test_R2)),
        col = "#B2182B20", border = NA)

legend("bottomright", legend = c("Training R²", "Test R² (±1 SD)"),
       col = c("#2166AC", "#B2182B"), pch = c(19, 17), lwd = 2, bty = "n")

# Add annotation
text(20, 0.6, "Gap = Overfitting\nConvergence = Sufficient data", cex = 0.9, adj = 0)

dev.off()
cat("Saved: output/day-4/26_learning_curve.png\n")

# -----------------------------------------------------------------------------
# COMPLETION
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("                  DAY 4 COMPLETE!\n")
cat("============================================================\n")
