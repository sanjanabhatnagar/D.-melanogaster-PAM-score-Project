# D.-melanogaster-PAM-score-Project

The script **`ExtractingcassettePSI_FromModENCODE.py`** extracts PSI values for cassette exons adjacent to introns listed in the PAM file, across modENCODE Drosophila tissues and developmental stages

Next, **`AddingSplicingdata_toPAMS.py`** script yields an output file with 
1. the PSI values for cassette exons adjacent to introns within a conserved RBP motif cluster identified in PAM score analysis.
2. Position of the intron relative to the candidate cassette exon - Upstream or Downstream of event.
3. Centered RBP expression across different tissue/developmental stages for which the PSI values have been reported. 
