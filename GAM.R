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


# Here -1 means there is no shared baseline PSI effect between events
# Here smooth functions were fitted separately for upstream and downstream intronic motif positions to capture nonlinear positional effects on splicing.
# Explicitly defined the link function as logit to keep PSI bound between 0 and 1, inspired from Zhao et al. 2013 (also, reflects true nature of PSI better
# The default link function is identity in gam! 

# In mgcv, random effects are implemented as penalized smooths. Hence, s(Meta_Corrected, bs = "re"), if I jus write Meta_Correctedf, it would be treated as a linear term 
# Check ?smooth.terms documentation in R 

fit_gamm_bypos <- gam(
  PSI ~ -1 + 
    s(RBP_centered, by = Upstream) +
    s(RBP_centered, by = Downstream) +
    s(Meta_Corrected, bs = "re"),
  family = quasibinomial(link = "logit"),
  data = data,
  method = "REML"
)

plot(fit_gamm_bypos, pages = 1, shade = TRUE)
summary(fit_gamm_bypos)
