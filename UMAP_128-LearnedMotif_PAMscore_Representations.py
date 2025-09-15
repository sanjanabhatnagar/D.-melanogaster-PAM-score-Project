import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns

pams_df = pd.read_csv('Finalrnarh_exon1_exon2_annotations_rnarh_interpretableN128L11_f1_b256_final_rep.csv', low_memory=False)

pams_df.drop(columns=['Unnamed: 0','FILE'], inplace=True)
pams_df_umap = pams_df.copy()

# The initial 7 columns contain annotations of introns. For instance, ordinal position of intron, type of intronic fragment - donor or acceptor, what kind of splicing event it flanks, etc.
# Hence, only processing the PAM score containing matrix.
pams_df_umap.iloc[:,7:] = pams_df_umap.iloc[:,7:].astype(float) 

# Centering the PAM scores on mean.
pams_df_umap.iloc[:,7:] = pams_df_umap.iloc[:,7:].apply(lambda x: x-x.mean())

# Generating UMAP 

# The columns containing intron annotations 
y_all = pams_df_umap[["meta_new", "ordinalposition", "meta_fragment", 'meta_type','gene','Exon_1', 'Exon_2']] 

# Just the PAM  score matrix, excluding all the annotation columns
X_all = pams_df_umap.drop(columns=["meta_new", "ordinalposition", "meta_fragment", 'meta_type','gene','Exon_1', 'Exon_2'])


reducer = umap.UMAP()
embedding_all = reducer.fit_transform(X_all)

umap_df_all = pd.DataFrame(embedding_all, columns=["UMAP1", "UMAP2"])

# Adding intron annotation data back to UMAP dataframe - Optional

umap_df_all["meta"] = y_all["meta_new"]
umap_df_all["Fragment"] = y_all["meta_fragment"]
umap_df_all["TypeofInt"] = y_all["meta_type"]
umap_df_all["ordinalposition"] = y_all["ordinalposition"]
umap_df_all["Exon_1"] = y_all["Exon_1"]
umap_df_all["Exon_2"] = y_all["Exon_2"]
umap_df_all["gene"] = y_all["gene"]

umap_df_all.to_csv('Finalrnarh_exon1_exon2_annotations_rnarh_interpretableN128L11_f1_b256_final_rep.csv') # This file is present in the associated release 1. Learned motif representations and dowsntream analysis data

# For generating the UMAP for instance, using Fragment annotations for the color of data points. 
plt.figure(figsize=(40, 20))

dark_palette = {
    " Acceptor 3'ss end " : "orange",
    " Donor 5'ss end ":"green",
    " Donor 5'ss and Acceptor 3'ss ":"blue"}


n_labels = umap_df_all["Fragment"].unique()

custom_palette = [dark_palette[label] for label in n_labels]  # order matters
palette = sns.color_palette(custom_palette[:len(n_labels)])


for i, label in enumerate(n_labels):
    subset = umap_df_all[umap_df_all["Fragment"] == label]
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

# Results in Figure 8A panel!
