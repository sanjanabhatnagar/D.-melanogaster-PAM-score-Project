import statsmodels.api as sm
import pandas as pd


X_with_const = sm.add_constant(X) # adds a column called const containing 1 in the matrix. Specifying the intercept if all idnependent variables are 0 or have no effect. So this is how we are setting the baseline.
# Here regressing UMAP1 dimension on all motif PAM scores.
model1 = sm.OLS(y_umap1, X_with_const).fit()
# Here regressing UMAP2 dimension on all motif PAM scores.
model2 = sm.OLS(y_umap2, X_with_const).fit()


results_umap1 = pd.DataFrame({
    "Motif": X_with_const.columns,
    "Coef_UMAP1": model1.params.values,
    "SE_UMAP1": model1.bse.values, # The standard error
    "t_UMAP1": model1.tvalues.values, # The t-stat value
    "pval_UMAP1": model1.pvalues.values
})

results_umap2 = pd.DataFrame({
    "Motif": X_with_const.columns,
    "Coef_UMAP2": model2.params.values,
    "SE_UMAP2": model2.bse.values,
    "t_UMAP2": model2.tvalues.values,
    "pval_UMAP2": model2.pvalues.values
})

coef_df = results_umap1.merge(results_umap2, on="Motif")
coef_df.to_csv('LinearRegression_RNARH_UMAP.csv')
