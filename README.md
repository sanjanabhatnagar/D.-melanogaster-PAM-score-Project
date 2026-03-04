# D.-melanogaster-PAM-score-Project

Prior to the analysis, we installed hal and ensured other conditions were met. Following repositories played a crucial role in our study.
 ```
gh repo clone ComparativeGenomicsToolkit/hal
gh repo clone pfenninglab/halLiftover-postprocessing
gh repo clone agduncan94/phylogenetic_augmentation_paper
 ```

## 1. Building D. melanogaster intronic sequences homolog sets by performing HAL liftover

We began by preparing D. melanogaster intronic sequences homolog sets from the other drosophilid species using custom scripts runLiftOverOnBed_gnu_parallel.sh which uses HAL liftover to fetch homologous sequences. We had also installed all the drosophilid whole genome fastas and drosophila_298.hal (Whole genome alignments). 

 ```

# bed file of intronic sequences of <100 bp length were prepared
# Strand-specific BED files processed separately to maintain directional integrity - Dmel6_60_intron_pv.bed and Dmel6_60_intron_nv.bed
# drosophila_species2.txt - contains names of species from the drosophilid family (homologous sequences are pulled from these species only)

# Running liftover for intronic sequences of genes on the positive strand
nohup bash /valr/sanjana/runLiftOverOnBed_gnu_parallel.sh /valr/sanjana/droso_pams_25/bed_files/Dmel6_60_intron_pv.bed D_MELANOGASTER ./droso_pams_25/drosophila_species2.txt dmel60_intronpv_ ./droso_pams_25/drosoInpv/ &
# Running liftover for intronic sequences of genes on the negative strand
nohup bash /valr/sanjana/runLiftOverOnBed_gnu_parallel.sh /valr/sanjana/droso_pams_25/bed_files/Dmel6_60_intron_nv.bed D_MELANOGASTER ./droso_pams_25/drosophila_species2.txt dmel60_intronnv_ ./droso_pams_25/drosoInnv/ &
 ```
This results in FASTA files segregated by species, each containing homologous sequences corresponding to the specified coordinates.

## 1.2. Multi-Species Homologous Sequence Grouping within the same fasta file 

The fasta files are accessed from drosoInpv/orthologs/per_species_fa/ and reverser complemented files (**`revcomp_pams.sh`**  which calls **`reverse_complement_pams.py`** ) from drosoInnv/orthologs/per_species_fa to the same folder to be processed further. 

 ```
nohup python convertFastaPerSpeciesToFastaPerSequence.py ./dm6_introns_bothStrands/ drosophila_chromosome_file.txt  &
 ```


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

## Re-running MAFFT on filtered sequence (with more than 70% of the sequences on the correct strand)

The files which did not pass the 0.7 threshold are removed after this step and the remaining files are realigned for downstream analysis.
Re-run MAFFT (but without adjusting directions this time) on the correct strandedness homolog sets obtained from the previous step to get accurately aligned homologs.

 ```
nohup bash run_multi_mafft_pams_strndcheck.sh drosophila_chromosome_file.txt 7 ./droso_introns_10exon_Feb2026/ > mafft_strndcheck.aln &
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

