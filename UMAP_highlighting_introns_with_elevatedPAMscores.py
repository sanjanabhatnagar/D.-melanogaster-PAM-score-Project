import matplotlib.pyplot as plt
import seaborn as sns

# Highlighting data points/introns with elevated PAM scores for strongest predictors of UMAP variation (From Linear Regression Results).

pams_df_matrix_1 = pams_df[(pams_df['8'].astype(float)>=3)]['meta_new'].copy() # matrix 8 
pams_df_matrix_2 = pams_df[(pams_df['22'].astype(float)>=3)]['meta_new'].copy() # matrix 22

mask_matrix =umap_df_all['meta'].isin(pams_df_matrix_1)
mask_matrix_2 =umap_df_all['meta'].isin(pams_df_matrix_2)

umap_df_all.loc[mask_matrix, 'matrix']='matrix 8'
umap_df_all.loc[mask_matrix_2, 'matrix']='matrix 22'

umap_df_all['matrix'].fillna('Others', inplace=True)


plt.figure(figsize=(40, 20))

dark_palette = {
    'matrix 8':'green',
    'matrix 22':'skyblue',
    'Others': 'lavender'
}


# Sort to make sure 'Others' is plotted first (i.e., in the back)
n_labels = ['Others'] + [label for label in umap_df_all["matrix"].unique() if label != 'Others']

custom_palette = [dark_palette[label] for label in n_labels]  # order matters
palette = sns.color_palette(custom_palette[:len(n_labels)])


#custom_palette = ['ghostwhite', 'red', 'ghostwhite', 'ghostwhite', 'orange']
#n_labels = umap_df["TypeofInt"].unique()
#palette = sns.color_palette(custom_palette[:len(n_labels)])
for i, label in enumerate(n_labels):
    subset = umap_df_all[umap_df_all["matrix"] == label]
    plt.scatter(
        subset["UMAP1"],
        subset["UMAP2"],
        edgecolor=palette[i],
        facecolor='none',
        linewidth=2,
        s=30,
        label=label
    )
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title("UMAP Projection of Motif Conservation Scores")
plt.xlim(-4, 20)
plt.xlabel('UMAP1', fontsize=24)
plt.ylabel('UMAP2', fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.grid(False)
plt.tight_layout()
plt.show()
