library(lme4)
library(lmerTest)
library(readr)
library(ggeffects)
library(ggplot2)

# Refer to http://r.qcbs.ca/workshop08/book-en/introduction-to-gams.html for more information on Generalized Additive Models.

df <- read_csv("/PAMS_final_Jan2026/1_Rbfox1_tissue_PSI.csv",
               show_col_types = FALSE) 

data <- subset(df, Position %in% c("Upstream", "Downstream"))

# This line conerts my meta_corrected/intron into a categorical variable.
data$Meta_Corrected <- factor(data$Meta_Corrected)

# fitting the data with a smoothed (non-linear) term 
# smooth terms are specified by expressions of the form s(x), where x is the non-linear predictor variable we want to smooth.
# Applied smooth function to interaction terms of (IR - Intronic Region; i is tissue or developmental stage context; j is given intronic fragment that belongs to the cluster with a conserved motif)
# (〖motif〗_(Upstream IR):RBP_expression )_ij  - FIXED EFFECT
# (〖motif〗_(Downstream IR):RBP_expression )_ij - FIXED EFFECT
# u_j+ ε_ij - RANDOM EFFECT

fit_gamm_bypos <- gam(
  PSI ~ -1 + 
    s(RBP_centered, by = Upstream, k = 10) +
    s(RBP_centered, by = Downstream, k = 10) +
    s(Meta_Corrected, bs = "re"),
  data = data,
  method = "REML"
)
plot(fit_gamm_bypos, pages = 1, shade = TRUE)
summary(fit_gamm_bypos)
