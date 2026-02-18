# D.-melanogaster-PAM-score-Project

The script **`ExtractingcassettePSI_FromModENCODE.py`** extracts PSI values for cassette exons adjacent to introns listed in the PAM file, across modENCODE Drosophila tissues and developmental stages

Next, **`AddingSplicingdata_toPAMS.py`** script yields an output file with 
1. the PSI values for cassette exons adjacent to introns within a conserved RBP motif cluster identified in PAM score analysis.
2. Position of the intron relative to the candidate cassette exon - Upstream or Downstream of event.
3. Centered RBP expression across different tissue/developmental stages for which the PSI values have been reported. 

```
Meta_Corrected	tissue	PSI	RBP	PAM	Position
7 | Acceptor 3'ss end | constitutive_downstream, skipped_upstream | FBtr0335228 | 8004864.0:8005009.0 | FBgn0037536 | CG2698	malphigian tubules	0.017026997	4.289800506	5.069569	Upstream
2 | Donor 5'ss end | constitutive_upstream(-), constitutive_totranscript_upstream(-), skipped_downstream(-) | FBtr0305011 | 17784096.0:17784537.0 | FBgn0262562 | CG43102	malphigian tubules	0.005725008	4.289800506	5.119941	Downstream
7 | Acceptor 3'ss end | constitutive_downstream, composite_alt_3_upstream | FBtr0076794 | 7535645.0:7535883.0 | FBgn0028582 | lqf	malphigian tubules	0.003253073	4.289800506	5.393003	Upstream
3 | Acceptor 3'ss end | constitutive_downstream(-), constitutive_totranscript_downstream(-), skipped_upstream(-) | FBtr0339079 | 15292970.0:15293426.0 | FBgn0265194 | Trpm	malphigian tubules	0.893981341	4.289800506	8.372929	Upstream
4 | Donor 5'ss end | constitutive_upstream, composite_alt_5_downstream | FBtr0305071 | 16292551.0:16292723.0 | FBgn0250823 | gish	malphigian tubules	0.006270104	4.289800506	6.003912	Downstream
2 | Donor 5'ss end | skipped_upstream(-), skipped_downstream(-) | FBtr0082206 | 10037891.0:10041971.0 | FBgn0053208 | Mical	malphigian tubules	0.885979932	4.289800506	3.926092	Downstream
4 | Donor 5'ss end | last_exon_upstream, skipped_downstream, skipped_upstream | FBtr0113007 | 21037995.0:21038745.0 | FBgn0031145 | Ntf-2	malphigian tubules	0.008290459	4.289800506	5.658145	Downstream
3 | Acceptor 3'ss end | constitutive_downstream(-), skipped_upstream(-) | FBtr0300713 | 1851497.0:1851661.0 | FBgn0035282 | CNMa	larval midgut	0.992498762	4.686368775	4.534145	Upstream
7 | Acceptor 3'ss end | constitutive_downstream, composite_alt_3_upstream | FBtr0076794 | 7535645.0:7535883.0 | FBgn0028582 | lqf	larval midgut	0.00842639	4.686368775	5.393003	Upstream
4 | Donor 5'ss end | constitutive_upstream, composite_alt_5_downstream | FBtr0305071 | 16292551.0:16292723.0 | FBgn0250823 | gish	larval midgut	0.002720449	4.686368775	6.003912	Downstream
2 | Donor 5'ss end | skipped_upstream(-), skipped_downstream(-) | FBtr0082206 | 10037891.0:10041971.0 | FBgn0053208 | Mical	larval midgut	0.796466417	4.686368775	3.926092	Downstream
4 | Donor 5'ss end | last_exon_upstream, skipped_downstream, skipped_upstream | FBtr0113007 | 21037995.0:21038745.0 | FBgn0031145 | Ntf-2	larval midgut	0.00188031	4.686368775	5.658145	Downstream
7 | Acceptor 3'ss end | constitutive_downstream, skipped_upstream | FBtr0335228 | 8004864.0:8005009.0 | FBgn0037536 | CG2698	larval fat	0.027468724	7.227830394	5.069569	Upstream
7 | Acceptor 3'ss end | constitutive_downstream, composite_alt_3_upstream | FBtr0076794 | 7535645.0:7535883.0 | FBgn0028582 | lqf	larval fat	0.015606387	7.227830394	5.393003	Upstream
3 | Acceptor 3'ss end | constitutive_downstream(-), constitutive_totranscript_downstream(-), skipped_upstream(-) | FBtr0339079 | 15292970.0:15293426.0 | FBgn0265194 | Trpm	larval fat	0.928365112	7.227830394	8.372929	Upstream
10 | Acceptor 3'ss end | skipped_upstream(-), skipped_downstream(-) | FBtr0342806 | 2499679.0:2500236.0 | FBgn0284408 | trol	larval fat	0.000778466	7.227830394	5.661129	Upstream
4 | Donor 5'ss end | constitutive_upstream, composite_alt_5_downstream | FBtr0305071 | 16292551.0:16292723.0 | FBgn0250823 | gish	larval fat	0.004866163	7.227830394	6.003912	Downstream
2 | Donor 5'ss end | skipped_upstream(-), skipped_downstream(-) | FBtr0082206 | 10037891.0:10041971.0 | FBgn0053208 | Mical	larval fat	0.840700038	7.227830394	3.926092	Downstream
......
```
