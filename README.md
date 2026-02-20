# D.-melanogaster-PAM-score-Project


## Running checks on MAFFT aligned files (--adjustdirectionaccurately parameter) to ensure all homologs are on same strand

The script **`strandedness_check_Feb2026.py`** computes the probabilities of 5'ss and 3'ss motifs being on sense strand vs. antisense strand for each homolog within a homolog set.

1. If more than 70% of the homologs have splice sites on the sense strand, the file is kept
2. If less than 70% of the homologs exhibit mixed strandedness and don't have the splice sites on the sense strand, the file is discarded.

Additionally this file also clips the first 10 bp exonic region for donor end fragments and last 10 bp exonic regions in D. melangaster and the corresponding aligned columns in all the other homologs. 


This script is run using the shell script **`STRND_CHECK_QC.sh`** 

The output files are - 
1. StrandednessCheck_QC_orthologgroups2L_aln.txt - Contains mean probabilities of splice sites being on sense vs. antisense strand across all homologs within a set.
2. Strn_check_fail_files_probablt70_2L_aln.txt - Contains names of files that did not pass 0.7 threshold for Mean(Motif|Sense) for both 5'ss and 3'ss motifs.

 StrandednessCheck_QC_orthologgroups2L_aln.txt - 
 ```
 FILE    MOTIF   Mean(Motif|Sense)       Stdev(Motif|Sense)      Total_sequences
2L_-_10001434_10001534.fa.aln   GTAAG   0.05224208861216475     0.10573835180575865     7
2L_-_10001434_10001534.fa.aln   CAG     0.699943776811742       0.2464845418349009      7
2L_-_10002087_10002187.fa.aln   GTAAG   0       0       176
2L_-_10002087_10002187.fa.aln   CAG     0.43004342156076025     0.34701081231300646     176	
......
```

## Adding ModEncode Analysis to PAM data

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

## Running Generalised Additive Mixed Model on PAM cassette events 

The **`GAMM.R`** script - models PSI using a generalized additive mixed model without a global intercept, incorporating separate nonlinear smooth effects of upstream and downstream intronic RBP signal and an event-specific random intercept to account for baseline splicing differences across events.

