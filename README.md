# Hammett Substituent Constant Statistical Analysis

## Project Overview

Hammett substituent constant statistical analysis of the para and meta positions for common EWG and EDGs in aromatic acids.

```
σ = log(K_X / K_H)
```



## Where Is the Conjugate Base?

A Hammett constant (σ) quantifies how much an aromatic substituent affects electron density, measuring its ability to donate or withdraw electrons (Evans, 2021).

- __σ > 0__: Electron-withdrawing group (EWG) - stabilizes negative charge
- __σ < 0__: Electron-donating group (EDG) - destabilizes negative charge
- __σ_meta vs σ_para__: Differences reveal resonance vs inductive effects

What happens when benzoic acid loses a proton?

```
        O               O⁻
        ‖               ‖
    ⬡—C—O—H    ⇌    ⬡—C—O⁻   +   H⁺
                         
  Benzoic Acid     Benzoate
   (neutral)      (negative)
```

The negative charge sits on the two oxygen atoms. Through resonance, this charge is delocalized across both oxygens. This is presented through the PI bond in the double bond.

How stable is this negative charge? The more stable the conjugate base, the more the equilibrium shifts toward deprotonation, and the stronger the acid.

## How Does Electron Withdrawal Stabilize Negative Charge?

Negative charge -> excess electron density -> repulsion, which is energetically unfavorable for a molecule. The molecule will attempt to distribute this charge.

Electron-withdrawing groups pull electron density away from the carboxylate, which leads to favorable charge distribution. This is the same reason why tertiary carbons are preferred in S<sub>N</sub>1.

```
         O⁻
         ‖
   EWG—⬡—C—O⁻
     ←←←

EWG pulls electron density 
away from the carboxylate,
spreading out the negative charge
```

Electron-donating groups do the opposite. They push more electron density toward the carboxylate, concentrating the negative charge further.

### Example

__Ethoxide (CH₃CH₂O⁻)__
- Alkyl group donates electron density.
- Charge stays localized on oxygen.
- Conjugate base is unstable ⇒ weak acid.

__Acetate (CH₃COO⁻)__
- Carbonyl withdraws electron density.
- Charge spreads over two oxygens by resonance.
- Conjugate base is stable ⇒ stronger acid.

Fundamentally it's because unstable conjugate bases do not like to give up their __positively__ charged proton, which sits right next to a highly __negative__ area.

## Why Meta vs Para Matters

```
           COOH
           |
           C 
        /  1  \
    6 C        C 2
      |        |
    5 C        C 3
        \  4  /
           C 
           
Position 2,6 = ortho
Position 5,3 = meta
Position 4 = para
```

There are __two different mechanisms__ for a substituent to affect the carboxylic acid.

### 1. Inductive Effect

```
            COOH
            |
            C ← ← ← 
          /   \      ↑
         C     C     |
         |     |     |  (electron pull
         C     C     |  through sigma bonds)
          \   /      |
           C ← ← ← ←
           |
          EWG
```

The inductive effect doesn't care about resonance structures, because it propogates through sigma bonds. It's simple to see the inductive effect should be stronger at the meta position because it is closer (has less bonds) to the carboxyl group.

### Resonance Effect

Resonance operates through the delocalized π electrons of the benzene ring. They only place charge at alternating locations, and for an EDG like NH<sub>2</sub>, the charge only occurs at the ortho and para positions, never at meta because you can push electrons to every other carbon (the exact method, I still need to research).

![wait](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEg6E_5PtQok5SYAe3-faWJjQWbwypJzkBK72LaoZwLMc22xuwwErUFPVfsXorcp55t3yXLil8PNv9rItcrFxaD7qG4LJ0vf1GVOWXSADZhBGrpft3RdjkfmzydlC53igrl24ID5RD3jtS_J/s1600/The-addition-of-an-electrophile-E-followed-by-elimination-of-a-proton-H-is-substitution-of-E-for-H-Electrophilic-Aromatic-Reactions-Reaction-of-Aromat.jpg "source: masterorganicchemistry.com")

The resonance donation creates negative charge density at the para carbon. 

The carboxylic acid is attached to a carbon that's attached to the para carbon. The electron density donated by NH₂ can increase electron density of the carboxylate group.

#### Note

Resonance puts negative charge on 2,4,6, and cannot place charge at the meta carbon. The meta position is "invisible" to resonance donation.

## Impact on Synthetic Pathways
- EWGs at benzylic position: stabilize carbocations → favor SN1
- EDGs through resonance: can stabilize both carbocations AND carbanions
- Understanding σ helps predict reaction selectivity

## Project Structure

```
hammett-analysis/
├── data/
│   ├── hammett_constants.csv     # Primary dataset (32 substituents)
│   ├── hammett_extended.csv      # Extended dataset with σ+ and σ-
│   └── hammett_clean.csv         # Processed dataset
├── scripts/
│   ├── 01_day1_setup.R           # Data loading and exploration
│   ├── 02_day2_eda.R             # Exploratory data analysis
│   ├── 03_day3_statistics.R      # Statistical tests
│   ├── 04_day4_modeling.R        # Predictive modeling
│   └── 05_day5_report.Rmd        # Final R Markdown report
├── output/
│   └── [generated plots and results]
├── docs/
│   └── [documentation]
└── README.md
```

## Data Sources

Primary data compiled from:
- C. Hansch and A. Leo, "Substituent Constants for Correlation Analysis in Chemistry and Biology," Wiley-Interscience, NY, 1979
- McDaniel & Brown, J. Org. Chem. 1958, 23, 420-427
- WiredChemist.com reference tables

## Variables in Dataset

| Variable | Description |
|----------|-------------|
| substituent | Chemical group formula |
| sigma_meta | Hammett constant (meta position) |
| sigma_para | Hammett constant (para position) |
| sigma_I | Inductive sigma constant |
| sigma_v | Charton's v (steric size) |
| pi | Hansch hydrophobicity parameter |
| Es | Taft steric parameter |
| MR | Molar refractivity |
| heteroatom | Key heteroatom (C, N, O, S, F, Cl, Br, I) |
| effect_type | EDG, EWG, neutral, weak_EDG, weak_EWG |
| functional_class | Chemical class (alkyl, halogen, amino, etc.) |

## Statistical Methods (Day 1-5)

1. __Exploratory Analysis__: Distributions, correlations, outlier detection
2. __Paired t-test__: Compare σ_meta vs σ_para
3. __ANOVA__: Test differences between heteroatom classes
4. __Linear Regression__: Predict σ_para from molecular descriptors
5. __Correlation Analysis__: Pearson/Spearman between sigma types

## Requirements

R >= 4.0.0 with packages:
- tidyverse
- ggplot2
- corrplot
- gridExtra
- knitr

## Running the Analysis

```r
# Set working directory to project root
setwd("path/to/hammett-analysis")

# Run Day 1 script
source("scripts/01_day1_setup.R")
```

## References

1. Hammett, L.P. (1937). "The Effect of Structure upon the Reactions of Organic Compounds." J. Am. Chem. Soc. 59(1): 96–103.
2. Hansch, C.; Leo, A.; Taft, R.W. (1991). "A Survey of Hammett Substituent Constants and Resonance and Field Parameters." Chem. Rev. 91(2): 165–195.