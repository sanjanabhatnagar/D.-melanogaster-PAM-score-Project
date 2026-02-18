import pandas as pd
modencode_expr_Mike = pd.read_csv('/logCPM_means.csv')

modencode_expr_Mike_RBP = modencode_expr_Mike[['Sample.Name','FBgn0052062']]

dev_stage= ['embryos, 0-2 hr','embryos,  2-4 hr','embryos, 4-6 hr','embryos, 6-8 hr','embryos, 8-10 hr','embryos, 10-12 hr','embryos, 12-14 hr','embryos, 14-16 hr','embryos, 16-18 hr','embryos, 18-20 hr','embryos, 20-22 hr','embryos, 22-24 hr','L1 stage larvae','L2 stage larvae','L3 stage larvae, 12 hr post-molt','WPP','pupae, 2 days after WPP','pupae, 3 days after WPP','pupae, 4 days after WPP','pupae, 12 hr after WPP','pupae, 24 hrs after WPP','adult female, 1 day after eclosion','adult female, 5 days after eclosion','adult female, 30 days after eclosion','adult male, 1 day after eclosion','adult male, 5 days after eclosion','adult male, 30 days after eclosion']
tissues = ['Adult Mated Female 1 day Post-eclosion Heads','Adult Mated Female 20 days Post-eclosion Heads','Adult Mated Female 4 days Post-eclosion Heads','Adult Mated Female 4 days Post-eclosion Ovaries','Adult Mated Male 1 day Post-eclosion Heads','Adult Mated Male 20 days Post-eclosion Heads','Adult Mated Male 4 days Post-eclosion AccessoryGlands','Adult Mated Male 4 days Post-eclosion Heads','Adult Mated Male 4 days Post-eclosion Testes','adult midgut','Adult Mixed Male and Female 1 day Post-eclosion Carcass','Adult Mixed Male and Female 1 day Post-eclosion Digestive System','Adult Mixed Male and Female 20 days Post-eclosion Carcass','Adult Mixed Male and Female 20 days Post-eclosion Digestive System','Adult Mixed Male and Female 4 days Post-eclosion Carcass','Adult Mixed Male and Female 4 days Post-eclosion Digestive System','Adult Virgin Female 1 day Post-eclosion Heads','Adult Virgin Female 20 days Post-eclosion Heads','Adult Virgin Female 4 days Post-eclosion Heads','Adult Virgin Female 4 days Post-eclosion Ovaries','L3 Carcass','L3 CNS','L3 Digestive System','L3 Fat body','L3 Imaginal Discs','L3 Salivary Glands','L3 stage larvae, clear gut PS(7-9)','L3 stage larvae, dark blue gut PS(1-2)','larval fat','larval midgut','malphigian tubules','WPP Fat Body','WPP Salivary Glands','WPP+2 days CNS']

modencode_expr_Mike_RBP_tissue = modencode_expr_Mike_RBP[modencode_expr_Mike_RBP['Sample.Name'].isin(tissues)]

modencode_expr_Mike_RBP_dev = modencode_expr_Mike_RBP[modencode_expr_Mike_RBP['Sample.Name'].isin(dev_stage)].copy()
modencode_expr_Mike_RBP_dev['Sample.Name'] = pd.Categorical(modencode_expr_Mike_RBP_dev['Sample.Name'],
                                                            categories=dev_stage,
                                                            ordered=True)
modencode_expr_Mike_RBP_dev = modencode_expr_Mike_RBP_dev.sort_values('Sample.Name')

modencode_expr_Mike_RBP_tissue.rename(columns={'FBgn0052062':'Rbfox1'}, inplace=True)
modencode_expr_Mike_RBP_dev.rename(columns={'FBgn0052062':'Rbfox1'}, inplace=True)


modencode_expr_Mike_RBP_tissue.to_csv('/modENCODE_analysis/Rbfox1_tissue_expression_logCPM.csv')
modencode_expr_Mike_RBP_dev.to_csv('/modENCODE_analysis/Rbfox1_dev_expression_logCPM.csv')
