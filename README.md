# D.-melanogaster-PAM-score-Project

The script **`ExtractingcassettePSI_FromModENCODE.py`** extracts PSI values for cassette exons adjacent to introns listed in the PAM file, across modENCODE Drosophila tissues and developmental stages

The script **`Processing_modEncodeExpressiondata.py`** extracts expression data for candidate RBP across modENCODE Drosophila tissues and developmental stages.

Next, **`AddingSplicingdata_toPAMS.py`** script yields an output file with 
1. the PSI values for cassette exons adjacent to introns within a conserved RBP motif cluster identified in PAM score analysis.
2. Position of the intron relative to the candidate cassette exon - Upstream or Downstream of event.
3. Centered RBP expression across different tissue/developmental stages for which the PSI values have been reported. 

```
                     Meta_Corrected	                                      tissue	          PSI	         RBP_centered	  PAM	   Position
(7 | Acceptor 3'ss end | constitutive_downstream, skipped_upstream |  malphigian tubules	0.017026997	   -2.929680759 5.069569	Upstream
8004864:8005009 | FBgn0037536 | CG2698)	
(2 | Donor 5'ss end | constitutive_upstream(-), skipped_downstream(-) | malphigian tubules	0.005725008	-2.929680759	5.119941	Downstream
17784096:17784537 | FBgn0262562 | CG43102)	
(3 | Acceptor 3'ss end | constitutive_downstream(-), skipped_upstream(-) | larval midgut	0.992498762	  -2.533112490	4.534145	Upstream
 1851497:1851661 | FBgn0035282 | CNMa)	
......
```

