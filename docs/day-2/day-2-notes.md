# Briefing

> Remember in Day 2, we are performing correlation comparisons, visualizing meta vs para, identifying outliers, and revisualizing in terms of functional group. Starting now, I will not be saying "inductive hammett constant" but rather sigma_I.

## Image 6, Correlation Heatmap

Certainly it is worth nothing the high-correlation top left corner. Particularly, the that sigma_I and sigma_meta have a 0.92 correlation, confirming sigma_meta measures primarily inductive effects.

You may nitpick the 0.08 not accounted for.
- Solvents (protic, aprotic, etc.)
- Electrostatic field interactions
- Secondary stereochemical interactions
- Experimental error (real possibility)

sigma_para and sigma_I have a 0.68 correlation, which is clearly weaker. We can confirm that in general, sigma_I does not cover resonance.

sigma_v and ES. These are two new parameters which are used in the Taft equation (a modification to the Hammett equation to account for spatial factors). In the molecules given, both of them measure molecular size.

They have a -0.96 correlation, which is strikingly close, and the negative sign is just a convention difference in how they're defined. 

![taft eqn](https://wikimedia.org/api/rest_v1/media/math/render/svg/1705ef410d2810a21e770b1ab2e034eb601d74f3 "source: Wikimedia Commons, 2008")

PI (Hansch hydrophobicity parameter) and sigma_meta have a tenuous downward correlation. EWGs tend to be polar and hydrophilic (negative PI), and EDGs like alkyl groups tend to be hydrophobic (positive PI). But EWGs and EDGs don't necessary have a strict mechanism of action, although most strong EDGs donate electrons primarily through mesomeric resonance.

### Core Logic
```R
# Calculate correlation matrix, handling missing values
cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

# Use layout() to create separate panels for heatmap and legend
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))

# image() plots the matrix as colored cells
image(1:ncol(cor_matrix), 1:nrow(cor_matrix), t(cor_matrix),
      col = color_palette, zlim = c(-1, 1))

# Add correlation values as text overlay
text(i, j, sprintf("%.2f", val), col = text_col)
```

## Image 7, Resonance Contribution Bar Chart

We straight up have an EDG that escaped the map (N(CH<sub>3</sub>)<sub>2</sub>). Jokes aside, we have each bar representing (sigma_para - sigma_meta), which isolates the resonance contribution.

Based on previous investigations, it's not hard to see the trend. Resonance donors (+R) are Lewis Bases that can donate a lone pair into an aromatic ring.

N(CH<sub>3</sub>)<sub>2</sub> dominates with -0.68, followed by NH<sub>2</sub> with -0.50. Both are known to be extremely strong nucleophiles. The amino groups are in general the strongest because nitrogen's lone pair is in a higher-energy orbital than oxygen's.

Positive values are resonance acceptors, and only NO shows strong resonance withdrawl (+0.29), and other EWGs like NO<sub>2</sub> and some carbonyls show slight positive values. Again, remember in Day 1 we showed resonance is a wild card in terms of strength, meta does have resonance.

The region close to the horizontal axis are primarily inductive alkyl groups, halogens, and first-degree halogenoalkanes. They lack strong resonance character, similar to meta. Although, fluorine (and to some extent, I and Br) break the trend.

Being the most electronegative element, F is a strong EWG, but in electrophilic aromatic substitution, halogens have considerable mesomeric donation effects, which partially counteracts inductive withdrawal. F is recognized as a weak ortho/para-director in certain reactions.

### Core Logic
```R
# Calculate resonance contribution
hammett$sigma_resonance <- hammett$sigma_para - hammett$sigma_meta

# Sort and assign colors based on threshold
hammett_sorted <- hammett[order(hammett$sigma_resonance), ]
bar_colors <- ifelse(hammett_sorted$sigma_resonance < -0.15, "#2166AC",  # donors
              ifelse(hammett_sorted$sigma_resonance > 0.15, "#B2182B",   # acceptors
                     "#999999"))                                         # neutral

# Create bar plot
barplot(hammett_sorted$sigma_resonance, names.arg = hammett_sorted$substituent,
        col = bar_colors, las = 2)
```

## Image 8, Resonance Contribution Bar Chart

Amino groups are the strongest EDG class having sigma_para = -0.75. We have seen this before, and it is primarily due to nitrogen's moderate electronegativity and strong mesomeric donation. Hydroxyl and alkoxy groups are weaker. sigma_para ≈ −0.26 to −0.37. Oxygen's higher electronegativity compared to nitrogen means stronger inductive withdrawal.

Halogens have a sigma_para = +0.18. Again, halogens are weak. Some notable EWGs are carbonyl, cyano, and nitro with a mean sigma_para = +0.46 to +0.78. Both inductive and resonance effects withdraw electrons.

Notice the gap between sigma_para and sigma_meta. The strength of EDGs comes primarily from resonance. In general, |sigma_para|>|sigma_meta| although ammonium are exceptions and groups in the middle like carbonyl groups are exceptions, where resonance withdrawal is present but less dramatic.

### Core Logic

```R
# Aggregate statistics by functional class
func_stats <- aggregate(cbind(sigma_meta, sigma_para) ~ functional_class, data = hammett, FUN = mean)

# Plot bars for sigma_meta
barplot(func_summary$mean_sigma_m, names.arg = classes, col = "#67A9CF")

# Overlay points for sigma_para
points(x_positions, func_summary$mean_sigma_p, pch = 18, col = "#B2182B")
```

## Image 9, Sigma-Plus vs Sigma-Minus

The sigma_+ and sigma_- are special Hammett constants that represent relative ability to stabilize charge in a reaction. sigma_+ measures positive charge (cations) and sigma_- measures negative charge (anions). The line is sigma_- = sigma_+.

Regular σ_para works for reactions where charge changes in the activated complex is negligible. sigma_+ is used for reactions that have significant charge changes in their steps, like S<sub>N</sub>1 and EA. Similarly, sigma_i is used for nucleophilic addition to carbonyls and some versions of E2.

### Key Observations
**1. EDGs appear in the bottom-left.**

N(CH₃)₂ has sigma_+ = −1.70, sigma_- = −0.83; NH₂ has sigma_+ = −1.30, sigma_- = −0.66. EDGs are excellent at stabilizing positive charge because they donate electrons to the empty p-orbital of carbocations via resonance. They are less effective at stabilizing negative charge, because on general sigma_- > sigma_+ for these groups, hence they are also to the left of the line.

**2. EWGs appear in the top-right.**

NO₂ has σ⁺ = +0.79, σ⁻ = +1.27; CN has sigma_+ = +0.66, sigma_- = +1.00. Now, the constants are positive, but they still work the same way. This measures the relative ability of EWGs to extract electrons. As you may expect these groups are very good at extracting negative charge, hence sigma_- > sigma_+ still holds.

Fluorine sits near the origin with sigma_+ = −0.07, sigma_- = +0.35. Inductive withdrawal and resonance donation nearly cancel for cations. Groups like fluorine (and to some extent other halogens) that are close to the origin may have some limited ability to stabilize charge. In fluorine, _inductive withdrawal wins for anions_.

### Core Logic

```R
# Plot with color coding by effect type
plot(hammett_ext_clean$sigma_plus, hammett_ext_clean$sigma_minus,
     col = ifelse(effect_type == "EDG", "#2166AC",
           ifelse(effect_type == "EWG", "#B2182B", "#666666")),
     pch = 19)

# Add diagonal reference line
abline(a = 0, b = 1, lty = 2)

# Label key substituents
text(sigma_plus[idx], sigma_minus[idx], labels = substituent, pos = 4)
```

## Image 10, Scatter Plot Matrix

Putting the Hammett constants together, we get a scatter plot matrix. The actual plots are on the bottom left, with their regression lines. The constants being compared are on the y=-x diagonal, and the Pearson Correlation Coefficients are on the top right.

In summary of our observations.
- sigma_para and sigma_meta have r = 0.9, which we already know from day 1 to be strong, but not perfect due to resonance.
- sigma_meta and sigma_I have a striking r = 0.92. We already concluded that the nearly perfect correlation confirms the meta position supports pure inductive effects. 
- sigma_para and sigma_I have r = 0.68. Inductive effects contribute, and the deviation comes from a word which occurs 25 times in the document, starts with r, ends with e.
- sigma_I and sigma_resonance have r = 0.15. They aren't really related, rather two different mechanisms of operation which may have trivial effects on each other.
    - A substituent can be inductive without being resonance-active (like CF₃), or  resonance-active without being inductive (eg. NH₂ vs OH).
- sigma_para and sigma_resonance have r = 0.82, which is expected, especially for strong EDGs whose main power comes from resonance.

### Strictly speaking, there are other factors aside from induction and resonance.
1. Steric inhibition of resonance (lone pair is not parallel to the ring)
2. Electrostatic field effects (E of a polar group can influence charges without passing through sigma bonds, influencing sigma_para)
3. Solvent effects (sometimes the ion-solvent complex can hinder effects)
4. Polarizability effects (easily polarizable LGs like iodine have easily distortable electron counds, hindering resonance donation)
5. Hyperconjugation (the methyl group have no lone pairs but C-H sigma may delocalize into the ring, this is NOT classical resonance but may affect sigma)
6. Measurement error (Hammett constants come from pKa, and this is open source data, which may not be 100% accurate)

### Core Logic
```R
# Custom panel functions for pairs plot
panel_scatter <- function(x, y) {
  points(x, y, pch = 19, col = "#2166AC80")
  abline(lm(y ~ x), col = "#B2182B", lwd = 2)
}

panel_cor <- function(x, y) {
  r <- cor(x, y, use = "pairwise.complete.obs")
  text(0.5, 0.5, sprintf("r = %.2f", r), cex = 1.5)
}

# Generate matrix of plots
pairs(plot_vars, lower.panel = panel_scatter, upper.panel = panel_cor)
```
