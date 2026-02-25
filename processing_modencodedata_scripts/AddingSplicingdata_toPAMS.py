import pandas as pd

df_cassette=pd.read_csv('/PAMS_final_Jan2026/PAMSfinaldm6_cassette.cdt')
RBP_cluster_df = pd.read_csv('/PAMS_final_Jan2026/PAMS_Motifs_processed/PAMS_Rbfox1_cluster.cdt', sep='\t')
RBP_expression_tissue =pd.read_csv('/modENCODE_analysis/Rbfox1_tissue_expression_logCPM.csv')

df_cassette_formerge = df_cassette.set_index('Meta_Corrected')
df_cassette_formerge_ = df_cassette_formerge.iloc[:,1:-8]
df_cassette_formerge_.drop(columns=['event_id'], inplace=True)
RBP_splicedata = pd.merge(RBP_cluster_df, df_cassette_formerge_, left_on='Meta_Corrected', right_index=True, how='left')
RBP_splicedata.set_index('Meta_Corrected', inplace=True)
df_splicing_data=RBP_splicedata.iloc[1:,-61:]

RBP_splicedata_tissue = df_splicing_data[['10_malphigian_tubules_median_psi', '11_larval_midgut_median_psi',
                                   '12_larval_fat_median_psi',
                                   '14_L3_stage_larvae_clear_gut_PS_7_9_median_psi','16_L3_Salivary_Glands_median_psi', '17_L3_Imaginal_Discs_median_psi',
                                   '18_L3_Fat_body_median_psi', '19_L3_Digestive_System_median_psi',
                                   '1_WPP_2_days_CNS_median_psi', '20_L3_CNS_median_psi',
                                   '21_L3_Carcass_median_psi','2_WPP_Salivary_Glands_median_psi','36_Adult_Virgin_Female_4_days_Post_eclosion_Ovaries_median_psi',
                                   '37_Adult_Virgin_Female_4_days_Post_eclosion_Heads_median_psi',
                                   '38_Adult_Virgin_Female_20_days_Post_eclosion_Heads_median_psi',
                                   '39_Adult_Virgin_Female_1_day_Post_eclosion_Heads_median_psi',
                                   '3_WPP_Fat_Body_median_psi',
                                   '40_Adult_Mixed_Male_and_Female_4_days_Post_eclosion_Digestive_System_median_psi',
                                   '41_Adult_Mixed_Male_and_Female_4_days_Post_eclosion_Carcass_median_psi',
                                   '42_Adult_Mixed_Male_and_Female_20_days_Post_eclosion_Digestive_System_median_psi',
                                   '43_Adult_Mixed_Male_and_Female_20_days_Post_eclosion_Carcass_median_psi',
                                   '44_Adult_Mixed_Male_and_Female_1_day_Post_eclosion_Digestive_System_median_psi',
                                   '45_Adult_Mixed_Male_and_Female_1_day_Post_eclosion_Carcass_median_psi',
                                   '46_adult_midgut_median_psi',
                                   '47_Adult_Mated_Male_4_days_Post_eclosion_Testes_median_psi',
                                   '48_Adult_Mated_Male_4_days_Post_eclosion_Heads_median_psi',
                                   '49_Adult_Mated_Male_4_days_Post_eclosion_AccessoryGlands_median_psi','50_Adult_Mated_Male_20_days_Post_eclosion_Heads_median_psi',
                                   '51_Adult_Mated_Male_1_day_Post_eclosion_Heads_median_psi',
                                   '52_Adult_Mated_Female_4_days_Post_eclosion_Ovaries_median_psi',
                                   '53_Adult_Mated_Female_4_days_Post_eclosion_Heads_median_psi',
                                   '54_Adult_Mated_Female_20_days_Post_eclosion_Heads_median_psi',
                                   '55_Adult_Mated_Female_1_day_Post_eclosion_Heads_median_psi','L13_L3_stage_larvae_dark_blue_gut_PS_1_2_median_psi']].copy()


RBP_splicedata_dev = df_splicing_data[['34_embryos_0_2_hr_median_psi',
                               '35_embryos__2_4_hr_median_psi','26_embryos_4_6_hr_median_psi','25_embryos_6_8_hr_median_psi','24_embryos_8_10_hr_median_psi','33_embryos_10_12_hr_median_psi','32_embryos_12_14_hr_median_psi',
                               '31_embryos_14_16_hr_median_psi','30_embryos_16_18_hr_median_psi','29_embryos_18_20_hr_median_psi','28_embryos_20_22_hr_median_psi','27_embryos_22_24_hr_median_psi','23_L1_stage_larvae_median_psi','22_L2_stage_larvae_median_psi','15_L3_stage_larvae_12_hr_post_molt_median_psi','4_WPP_median_psi',
                               '8_pupae_2_days_after_WPP_median_psi','6_pupae_3_days_after_WPP_median_psi','5_pupae_4_days_after_WPP_median_psi',
                               '9_pupae_12_hr_after_WPP_median_psi','7_pupae_24_hrs_after_WPP_median_psi',
                               '58_adult_male_1_day_after_eclosion_median_psi','61_adult_female_1_day_after_eclosion_median_psi',
                               '56_adult_male_5_days_after_eclosion_median_psi','59_adult_female_5_days_after_eclosion_median_psi',
                               '57_adult_male_30_days_after_eclosion_median_psi','60_adult_female_30_days_after_eclosion_median_psi']].copy()

RBP_splicedata_dev.rename(columns={'34_embryos_0_2_hr_median_psi':'embryos, 0-2 hr', '35_embryos__2_4_hr_median_psi':'embryos,  2-4 hr', '26_embryos_4_6_hr_median_psi':'embryos, 4-6 hr', '25_embryos_6_8_hr_median_psi':'embryos, 6-8 hr', '24_embryos_8_10_hr_median_psi':'embryos, 8-10 hr', '33_embryos_10_12_hr_median_psi':'embryos, 10-12 hr', '32_embryos_12_14_hr_median_psi':'embryos, 12-14 hr', '31_embryos_14_16_hr_median_psi':'embryos, 14-16 hr', '30_embryos_16_18_hr_median_psi':'embryos, 16-18 hr', '29_embryos_18_20_hr_median_psi':'embryos, 18-20 hr', '28_embryos_20_22_hr_median_psi':'embryos, 20-22 hr', '27_embryos_22_24_hr_median_psi':'embryos, 22-24 hr', '23_L1_stage_larvae_median_psi':'L1 stage larvae', '22_L2_stage_larvae_median_psi':'L2 stage larvae', '15_L3_stage_larvae_12_hr_post_molt_median_psi':'L3 stage larvae, 12 hr post-molt', '4_WPP_median_psi':'WPP', '8_pupae_2_days_after_WPP_median_psi':'pupae, 2 days after WPP', '6_pupae_3_days_after_WPP_median_psi':'pupae, 3 days after WPP', '5_pupae_4_days_after_WPP_median_psi':'pupae, 4 days after WPP', '9_pupae_12_hr_after_WPP_median_psi':'pupae, 12 hr after WPP', '7_pupae_24_hrs_after_WPP_median_psi':'pupae, 24 hrs after WPP', '58_adult_male_1_day_after_eclosion_median_psi':'adult male, 1 day after eclosion', '61_adult_female_1_day_after_eclosion_median_psi':'adult female, 1 day after eclosion', '56_adult_male_5_days_after_eclosion_median_psi':'adult male, 5 days after eclosion', '59_adult_female_5_days_after_eclosion_median_psi':'adult female, 5 days after eclosion', '57_adult_male_30_days_after_eclosion_median_psi':'adult male, 30 days after eclosion', '60_adult_female_30_days_after_eclosion_median_psi':'adult female, 30 days after eclosion'}, inplace=True)
RBP_splicedata_tissue.rename(columns={'10_malphigian_tubules_median_psi':'malphigian tubules', '11_larval_midgut_median_psi':'larval midgut', '12_larval_fat_median_psi':'larval fat', '14_L3_stage_larvae_clear_gut_PS_7_9_median_psi':'L3 stage larvae, clear gut PS(7-9)', '16_L3_Salivary_Glands_median_psi':'L3 Salivary Glands', '17_L3_Imaginal_Discs_median_psi':'L3 Imaginal Discs', '18_L3_Fat_body_median_psi':'L3 Fat body', '19_L3_Digestive_System_median_psi':'L3 Digestive System', '1_WPP_2_days_CNS_median_psi':'WPP+2 days CNS', '20_L3_CNS_median_psi':'L3 CNS', '21_L3_Carcass_median_psi':'L3 Carcass', '2_WPP_Salivary_Glands_median_psi':'WPP Salivary Glands', '36_Adult_Virgin_Female_4_days_Post_eclosion_Ovaries_median_psi':'Adult Virgin Female 4 days Post-eclosion Ovaries', '37_Adult_Virgin_Female_4_days_Post_eclosion_Heads_median_psi':'Adult Virgin Female 4 days Post-eclosion Heads', '38_Adult_Virgin_Female_20_days_Post_eclosion_Heads_median_psi':'Adult Virgin Female 20 days Post-eclosion Heads', '39_Adult_Virgin_Female_1_day_Post_eclosion_Heads_median_psi':'Adult Virgin Female 1 day Post-eclosion Heads', '3_WPP_Fat_Body_median_psi':'WPP Fat Body', '40_Adult_Mixed_Male_and_Female_4_days_Post_eclosion_Digestive_System_median_psi':'Adult Mixed Male and Female 4 days Post-eclosion Digestive System', '41_Adult_Mixed_Male_and_Female_4_days_Post_eclosion_Carcass_median_psi':'Adult Mixed Male and Female 4 days Post-eclosion Carcass', '42_Adult_Mixed_Male_and_Female_20_days_Post_eclosion_Digestive_System_median_psi':'Adult Mixed Male and Female 20 days Post-eclosion Digestive System', '43_Adult_Mixed_Male_and_Female_20_days_Post_eclosion_Carcass_median_psi':'Adult Mixed Male and Female 20 days Post-eclosion Carcass', '44_Adult_Mixed_Male_and_Female_1_day_Post_eclosion_Digestive_System_median_psi':'Adult Mixed Male and Female 1 day Post-eclosion Digestive System', '45_Adult_Mixed_Male_and_Female_1_day_Post_eclosion_Carcass_median_psi':'Adult Mixed Male and Female 1 day Post-eclosion Carcass', '46_adult_midgut_median_psi':'adult midgut', '47_Adult_Mated_Male_4_days_Post_eclosion_Testes_median_psi':'Adult Mated Male 4 days Post-eclosion Testes', '48_Adult_Mated_Male_4_days_Post_eclosion_Heads_median_psi':'Adult Mated Male 4 days Post-eclosion Heads', '49_Adult_Mated_Male_4_days_Post_eclosion_AccessoryGlands_median_psi':'Adult Mated Male 4 days Post-eclosion AccessoryGlands', '50_Adult_Mated_Male_20_days_Post_eclosion_Heads_median_psi':'Adult Mated Male 20 days Post-eclosion Heads', '51_Adult_Mated_Male_1_day_Post_eclosion_Heads_median_psi':'Adult Mated Male 1 day Post-eclosion Heads', '52_Adult_Mated_Female_4_days_Post_eclosion_Ovaries_median_psi':'Adult Mated Female 4 days Post-eclosion Ovaries', '53_Adult_Mated_Female_4_days_Post_eclosion_Heads_median_psi':'Adult Mated Female 4 days Post-eclosion Heads', '54_Adult_Mated_Female_20_days_Post_eclosion_Heads_median_psi':'Adult Mated Female 20 days Post-eclosion Heads', '55_Adult_Mated_Female_1_day_Post_eclosion_Heads_median_psi':'Adult Mated Female 1 day Post-eclosion Heads', 'L13_L3_stage_larvae_dark_blue_gut_PS_1_2_median_psi':'L3 stage larvae, dark blue gut PS(1-2)'}, inplace=True)

RBP_splicedata_tissue_reset = RBP_splicedata_tissue.reset_index()
RBP_splicedata_tissue_long = RBP_splicedata_tissue_reset.melt(id_vars='Meta_Corrected',
                                                              var_name='tissue',
                                                              value_name='PSI')


RBP_expression_tissue.drop(columns = 'Unnamed: 0', inplace=True)
Perintron_splicedata_pertissue_ = pd.merge(RBP_splicedata_tissue_long, RBP_expression_tissue, left_on='tissue', right_on='Sample.Name', how='left')
Perintron_splicedata_pertissue = pd.merge(Perintron_splicedata_pertissue_, RBP_cluster_df[['Meta_Corrected','Rbfox1']], left_on='Meta_Corrected', right_on='Meta_Corrected', how='left')
Perintron_splicedata_pertissue.drop(columns='Sample.Name', inplace=True)
Perintron_splicedata_pertissue.dropna(subset=['PSI'], axis=0, inplace=True)
Perintron_splicedata_pertissue.columns = ['Meta_Corrected', 'tissue', 'PSI','RBP', 'PAM']
mask0 = (Perintron_splicedata_pertissue['Meta_Corrected'].astype(str).str.contains('skipped_upstream|composite_alt_5_upstream|composite_alt_3_upstream') & Perintron_splicedata_pertissue['Meta_Corrected'].astype(str).str.contains('Acceptor'))

mask1 = (Perintron_splicedata_pertissue['Meta_Corrected'].astype(str).str.contains('skipped_downstream|composite_alt_5_downstream|composite_alt_3_downstream') & Perintron_splicedata_pertissue['Meta_Corrected'].astype(str).str.contains('Donor') )



Perintron_splicedata_pertissue.loc[mask0, 'Position']= 'Upstream'
Perintron_splicedata_pertissue.loc[mask1, 'Position']= 'Downstream'
Perintron_splicedata_pertissue['PAM_centered'] = Perintron_splicedata_pertissue['PAM'] - Perintron_splicedata_pertissue['PAM'].mean()
Perintron_splicedata_pertissue['RBP_centered'] = Perintron_splicedata_pertissue['RBP'] - Perintron_splicedata_pertissue['RBP'].mean()

Perintron_splicedata_pertissue['PAM*RBP_tissue'] = Perintron_splicedata_pertissue['PAM_centered'] * Perintron_splicedata_pertissue['RBP_centered']

Perintron_splicedata_pertissue['PAM_RBP']=Perintron_splicedata_pertissue['PAM'] * Perintron_splicedata_pertissue['RBP']
Perintron_splicedata_pertissue['group'] = 'conserved'

Perintron_splicedata_pertissue.dropna(subset='Position', inplace=True)
Perintron_splicedata_pertissue.to_csv('/PAMS_final_Jan2026/1_Rbfox1_tissue_PSI.csv')
