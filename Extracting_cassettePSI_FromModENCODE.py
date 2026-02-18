import pandas as pd

pams_df = pd.read_csv('/PAMS_final_Jan2026/PAMS_final_dm6_introns.cdt', sep='\t', low_memory=False)
pams_df = pams_df.iloc[2:,:]
pams_df['Gene_ID'] = pams_df['Meta_Corrected'].astype(str).str.split('|').str[5]
pams_df['start'] = pams_df['NAME'].astype(str).str.split('_').str[2]
pams_df['end'] = pams_df['NAME'].astype(str).str.split('_').str[3]
pams_df['end']=pams_df['end'].astype(str).str.replace('.fa','')
pams_df[['Chr','strand']]=pams_df['NAME'].astype(str).str.split('_',expand=True)[[0,1]]
pams_df[['start','end']] = pams_df[['start','end']].astype(int)

pams_df['Gene_ID'] = pams_df['Gene_ID'].astype(str).str.replace(' ','')
pams_df['ID'] = pams_df['NAME'].astype(str).str.split('.fa').str[0]

heatmap_unfiltered = pd.read_csv('/Modencode_processed_splicingdata_Mike/cassette.tsv',sep='\t')
heatmap_unfiltered['gene_id']= heatmap_unfiltered['gene_id'].astype(str).str.replace('gene:','')
heatmap_unfiltered[['exon_start', 'exon_end']] = heatmap_unfiltered['spliced_with_coord'].astype(str).str.split('-', expand=True)
heatmap_unfiltered[['exon_start', 'exon_end']]= heatmap_unfiltered[['exon_start', 'exon_end']].astype(int)

def select_low_variance_event(group):
    variance_cols = [col for col in group.columns if 'var_psi' in col]
    if group['event_id'].nunique() == 1:
        return group
    else:
        event_variances = group.groupby('event_id')[variance_cols].mean().mean(axis=1) # mean of mean per sample variance
        best_event_id = event_variances.idxmin()
        return group[group['event_id'] == best_event_id]

heatmap = heatmap_unfiltered.groupby(['exon_start', 'exon_end'], group_keys=False).apply(select_low_variance_event)
heatmap.drop(columns=['exon_start', 'exon_end'], inplace=True)
heatmap = heatmap.drop(
    columns=[col for col in heatmap.columns if 'var_psi' in col]
)

def cassettePSI_from_lsv_junc(group):
    group.iloc[:,18:] = group.iloc[:,18:].replace(0, 1e-9)
    alt_exon = group[((group['junction_name'] == 'C1_A') | (group['junction_name'] == 'C2_A'))]
    PSI = alt_exon.iloc[:,18:].mean()
    return PSI

heatmap_cassette_PSI = heatmap.groupby('event_id').apply(cassettePSI_from_lsv_junc)
heatmap_cassette_PSI.reset_index(inplace=True)
casette_PSI=pd.merge(heatmap_cassette_PSI, heatmap.iloc[:,1:14], on='event_id', how='left')
casette_PSI[['exon_start', 'exon_end']] = casette_PSI['spliced_with_coord'].astype(str).str.split('-', expand=True)
casette_PSI[['exon_start', 'exon_end']]= casette_PSI[['exon_start', 'exon_end']].astype(int)
casette_PSI = casette_PSI[(casette_PSI['spliced_with']=='A')]
psi_cols = [col for col in casette_PSI.columns if col.endswith('_median_psi')]
casette_PSI[psi_cols] = casette_PSI[psi_cols].apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0)

def splice_data_map(row):
    mod_row = None
    pam_strand = row['strand']
    pam_meta = row['Meta_Corrected']
    pam_start = row['start']
    pam_end = row['end']
    pam_gene = row['Gene_ID']
    if pam_strand == '+':
        if pd.notna(pam_meta) and 'Donor 5\'ss and Acceptor 3\'ss' in pam_meta:#in pam_meta
            match = casette_PSI[(casette_PSI['exon_end'] == pam_start - 1) & (casette_PSI['gene_id'] == pam_gene) & (casette_PSI['exon_start'] == pam_end + 1)]
        else:
            match = casette_PSI[((casette_PSI['exon_end'] == pam_start - 1) & (casette_PSI['gene_id'] == pam_gene)) | ((casette_PSI['exon_start'] == pam_end + 1) & (casette_PSI['gene_id'] == pam_gene))]

        if not match.empty:
            mod_row = match.copy()
            mod_row[['ID','strand','start','end','Meta_Corrected']] = row[['ID','strand','start','end','Meta_Corrected']].values


    elif pam_strand == '-':
        if pd.notna(pam_meta) and 'Donor 5\'ss and Acceptor 3\'ss' in pam_meta:          #'Donor 5\'ss and Acceptor 3\'ss' in pam_meta
            match =casette_PSI[(casette_PSI['exon_end'] == pam_start - 1) & (casette_PSI['gene_id'] == pam_gene) & (casette_PSI['exon_start'] == pam_end + 1)]
        else:
            match = casette_PSI[((casette_PSI['exon_end'] == pam_start - 1) & (casette_PSI['gene_id'] == pam_gene)) | ((casette_PSI['exon_start'] == pam_end + 1) & (casette_PSI['gene_id'] == pam_gene))]

        if not match.empty:
            mod_row = match.copy()
            mod_row[['ID','strand','start','end','Meta_Corrected']] = row[['ID','strand','start','end','Meta_Corrected']].values

    return mod_row

df_splice = pams_df.apply(splice_data_map, axis=1)
df_splice = pd.concat([df for df in df_splice], ignore_index=True)
df_splice = df_splice.drop(
    columns=[col for col in df_splice.columns if 'var_psi' in col]
)
df_splice_cleaned_1 = df_splice.drop(columns=['complex','denovo','reference_exon_coord','seqid', 'strand', 'lsv_id','spliced_with_coord','junction_name','junction_coord']).drop_duplicates()
df_cassette = df_splice_cleaned_1

df_cassette.to_csv('/PAMS_final_Jan2026/PAMSfinaldm6_cassette.cdt')
