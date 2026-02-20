import sys
import os
import pandas as pd
import numpy as np
import math
from Bio.Seq import Seq
from Bio import motifs
from itertools import groupby

input_aln_fasta = sys.argv[1]
output_file = sys.argv[2]
output_file_probablt70 = sys.argv[3]
type_=sys.argv[4]
option_3=None
adjusted_coords_file = './Dmel6_60_intron_both_+10exons.bed.csv'
name, ext = os.path.splitext(input_aln_fasta)
input_aln_fasta_modified = f"{name}_exons_clipped{ext}"
to_be_discarded=[]
# The following code removes the dashes from alignment and ensures each sequence is present in a single line
#Basically it creates a dictionary with the header >D_MEL etc. and the following sequence as the value!

def clip_from_start(seq, len_2_clip):
    count_nts = 0
    i = 0
    while i < len(seq) and count_nts < len_2_clip:
        if seq[i] != '-':
            count_nts += 1
        i += 1
    return i


def clip_from_end(seq, len_2_clip):
    count_nts = 0
    i = len(seq) - 1
    while i >= 0 and count_nts < len_2_clip:
        if seq[i] != '-':
            count_nts += 1
        i -= 1
    return i

def clip_exons(alignment, coord_df):
    clipped_alignment = {}
    #>D_MELANOGASTER::4:261128-261239
    # part 0 will be D_MELANOGASTER, part 1 will be chromosome, part 2 will be start-end coords
    dmel_header = [h for h in alignment if h.startswith('>D_MELANOGASTER')]
    parts = dmel_header[0].split(':')
    chrom = parts[2]
    start, end = map(int, parts[3].split('-'))
    # Then I try to find where these chr, start and end match in the coord_Df
    coord_df['start'] = coord_df['start'].astype(int)
    coord_df['chr'] = coord_df['chr'].astype(str)
    coord_df['end'] = coord_df['end'].astype(int)
    match = coord_df[(coord_df['chr'] == chrom) & (coord_df['start'] == start) & (coord_df['end']== end)]
    if match.empty:#("Highly unlikely but still")
        print("Error: No matching coordinates found in the coordinate file!")
        return
    typeofend = match.iloc[0]['Type']
    strand = match.iloc[0]['strand']
    parts_=match.iloc[0]['ID'].split('_')
    start_real = parts_[2]
    end_real=parts_[3]

    # So now, code first fetches 10 bp and the corresponding alignment columns from the alignment for the dmel sequence only (our ref)
    dmel_header = [h for h in alignment if h.startswith('>D_MELANOGASTER')][0] # Take the header
    dmel_seq = alignment[dmel_header] # Find the seq for the header

    if typeofend == "Acceptor 3'ss end":
        aln_col_position_accep = clip_from_end(dmel_seq, 10)

    elif typeofend == "Donor 5'ss end":
        aln_col_position_donor = clip_from_start(dmel_seq, 10)

    elif typeofend == "Donor 5'ss and Acceptor 3'ss":
        aln_col_position_donor = clip_from_start(dmel_seq, 10)
        aln_col_position_accep = clip_from_end(dmel_seq, 10)

    # Once the columns are fetched, now I can clip all seqs at once whether they are ref or other species
    for header, sequence in alignment.items():
        if typeofend == "Acceptor 3'ss end":
            clipped_seq = sequence[:aln_col_position_accep+1]
            print(f"{header}\n{clipped_seq}")

        elif typeofend == "Donor 5'ss end":
            clipped_seq = sequence[aln_col_position_donor:]
            print(f"{header}\n{clipped_seq}")

        elif typeofend == "Donor 5'ss and Acceptor 3'ss":
            clipped_seq = sequence[aln_col_position_donor:aln_col_position_accep+1]
            print(f"{header}\n{clipped_seq}")

        ## Finally I also modify the headers so that I can distinguish sequences with clipped files from input files!!
        if header.startswith('>D_MELANOGASTER'):
            clipped_header = f'>D_MELANOGASTER::{chrom}:{start_real}-{end_real}'
        else:
            clipped_header = f'{header}::{typeofend}:coords_not_updated'

        clipped_alignment[clipped_header] = clipped_seq

    with open(input_aln_fasta_modified, 'w') as temp:
        for header, sequence in clipped_alignment.items():
            temp.write(f'{header}\n{sequence}\n')


def complement(s):
    comp_dict = {'A': 'T', 'a': 't', 'G': 'C', 'g': 'c', 'C': 'G', 'c': 'g', 'T': 'A', 't': 'a'}
    return ''.join(comp_dict.get(base, '') for base in s)
def reverse_complement(func, z):
    return func(z)[::-1]

#These PFMs are from Mount et al. 1992, drosophila splice site paper!

# G T R A G T Drosophila (all introns) 5' splice site sequences Table 2A, Mount et al. 1992
pfm_5 = [
    [0, 0, 60, 71, 9, 11],
    [0, 0, 1, 9, 2, 14],
    [100, 0, 35, 9, 82, 6],
    [0, 100, 4, 11, 6, 68]
]
pfm_dict = {
    "A": pfm_5[0],
    "C": pfm_5[1],
    "G": pfm_5[2],
    "T": pfm_5[3]
}

PFM_5ss = motifs.Motif(counts=pfm_dict)

# Y T T _ C A G Drosophila (all introns) 3' splice site sequences Table 2B, Mount et al. 1992
pfm_3 = [
    [19, 10, 11, 28, 5, 100, 0],
    [36, 28, 20, 23, 68, 0, 0],
    [6, 5, 4, 23, 0, 0, 100],
    [40, 57, 64, 25, 27, 0, 0]
]
pfm3_dict = {
    "A": pfm_3[0],
    "C": pfm_3[1],
    "G": pfm_3[2],
    "T": pfm_3[3]
}
PFM_3ss = motifs.Motif(counts=pfm3_dict)
pwm5 = PFM_5ss.counts.normalize(pseudocounts=0.1)
pwm_5 = pwm5.log_odds()
pwm3 = PFM_3ss.counts.normalize(pseudocounts=0.1)
pwm_3 = pwm3.log_odds()

def scan(pssm, seq, minscore, motif_id):
    results = []
    for position, score in pssm.search(seq, threshold=minscore, both=False):
        end_position = position + len(pssm.consensus)
        values = [motif_id,
                  position + 1, end_position,
                  str(seq[position:end_position]),  # .transcribe()
                  round(score, 3)]
        results.append(values)

    if not results:
        return float('-inf')

    best_match = max(
        results,
        key=lambda x: x[-1] if not math.isnan(x[-1]) else float('-inf')
    ) # Sometimes multiple sites are matched to the splice site because I have a lower threshold to get the scores of splice sites on antisense/reverse complement strands.
    score = best_match[4]
    site = best_match[3]
    return score

# I hope this is the right way to convert log-odds score to probabilities !
def sigmoid(score):
    return 1 / (1 + np.exp(-score))

def strand_check_(aln):
    c = 0
    motif5_sense_P_seqs =[]
    motif3_sense_P_seqs = []

    for seq in aln.values():
        seq1 = seq.upper()
        seq1 =seq.strip()
        seq1 = seq.replace('-','')
        if not (seq.startswith('>')):
            c+=1
            seq_ = Seq(seq1)
            seq_rc = reverse_complement(complement, seq_)

            if len(seq_) >= 16: # That's what I am specifying because the sequence should be larger than the splice site PWM lengths!!!
                # Here I am getting log-odds score for 5'ss and 3'ss in sequence
                motif5ss_seq_score= scan(pwm_5, seq_[:7], -20,'5ss') # G T R A G T Drosophila (all introns) 5' splice site sequences Table 2A, Mount et al. 1992
                motif3ss_seq_score = scan(pwm_3, seq_[-10:], -20,'3ss') #  Y T T _ C A G Drosophila (all introns) 3' splice site sequences Table 2B, Mount et al. 1992

                # Here I am getting log-odds score for 5'ss and 3'ss in reverse_complement of the sequence
                motif5ss_seqrc_score = scan(pwm_5, seq_rc[:7], -20,'5ss') # G T R A G T Drosophila (all introns) 5' splice site sequences Table 2A, Mount et al. 1992
                motif3ss_seqrc_score = scan(pwm_3, seq_rc[-10:], -20,'3ss') #  Y T T _ C A G Drosophila (all introns) 3' splice site sequences Table 2B, Mount et al. 1992

            else:
                motif5ss_seq_score = 0
                motif3ss_seq_score = 0
                motif5ss_seqrc_score = 0
                motif3ss_seqrc_score = 0


            #Converting the scores to probabilities using sigmoid function
            motif5ss_seq_P = sigmoid(motif5ss_seq_score)
            motif3ss_seq_P = sigmoid(motif3ss_seq_score)

            motif5ss_seqrc_P = sigmoid(motif5ss_seqrc_score)
            motif3ss_seqrc_P = sigmoid(motif3ss_seqrc_score)

            #Calculating the probability of finding 5'ss or 3'ss motif given the strand is +ve/sense
            motif5_sense_P = motif5ss_seq_P / (motif5ss_seq_P + motif5ss_seqrc_P)
            motif3_sense_P = motif3ss_seq_P / (motif3ss_seq_P + motif3ss_seqrc_P)

            motif5_sense_P_seqs.append(motif5_sense_P)
            motif3_sense_P_seqs.append(motif3_sense_P)


    mean_5ss_score = 0 if np.isnan(np.mean(motif5_sense_P_seqs)) else np.mean(motif5_sense_P_seqs)
    std_5ss_score  = 0 if np.isnan(np.std(motif5_sense_P_seqs))  else np.std(motif5_sense_P_seqs)

    mean_3ss_score = 0 if np.isnan(np.mean(motif3_sense_P_seqs)) else np.mean(motif3_sense_P_seqs)
    std_3ss_score  = 0 if np.isnan(np.std(motif3_sense_P_seqs))  else np.std(motif3_sense_P_seqs)

    return(mean_5ss_score, std_5ss_score,mean_3ss_score, std_3ss_score, c)

with open(input_aln_fasta, 'r') as handle:
    aln_ = handle.readlines()

alignments= {}
header = None
for is_header, group in groupby(aln_, key=lambda x: x.startswith('>')):
    if is_header:
        header = next(group).strip()
    else:
        sequence = ''.join(line.strip() for line in group)
        alignments[header] = sequence

if type_=='introns':
    if option_3=='ref_only':
        dmel_only = {key: alignments[key] for key in alignments if "D_MELANOGASTER" in key}
        print("Working on D_MEL")
        avg_5motif_sense, std_5motif_sense, avg_3motif_sense, std_3motif_sense, total_seq = strand_check_(dmel_only)
    elif option_3==None:
        avg_5motif_sense, std_5motif_sense,avg_3motif_sense, std_3motif_sense, total_seq = strand_check_(alignments)

elif type_=='exons':
    cdf = pd.read_csv(adjusted_coords_file,sep='\t')
    clip_exons(alignments, cdf)

    with open(input_aln_fasta_modified, 'r') as handle:
        aligns_2 = handle.readlines()
    alignments_2={}

    for line in aligns_2:
        if line.startswith('>'):
            header = line.strip()
        elif line:
            sequence = line.strip()
            alignments_2[header] = sequence

    if option_3=='ref_only':
        dmel_only = {key: alignments_2[key] for key in alignments_2 if "D_MELANOGASTER" in key}
        print("Working on D_MEL")
        print("in here")
        print(alignments_2)
        avg_5motif_sense, std_5motif_sense, avg_3motif_sense, std_3motif_sense, total_seq = strand_check_(dmel_only)

    elif option_3==None:
        avg_5motif_sense, std_5motif_sense,avg_3motif_sense, std_3motif_sense, total_seq = strand_check_(alignments_2)
        if (avg_5motif_sense<0.7 and avg_3motif_sense<0.7):
            to_be_discarded.append(f'{os.path.basename(input_aln_fasta_modified)}')


with open(output_file, 'a') as handle2:
    handle2.write(f'{os.path.basename(input_aln_fasta)}\tGTAAG\t{avg_5motif_sense}\t{std_5motif_sense}\t{total_seq}\n')
    handle2.write(f'{os.path.basename(input_aln_fasta)}\tCAG\t{avg_3motif_sense}\t{std_3motif_sense}\t{total_seq}\n')

with open(output_file_probablt70, 'a') as handle3:
    for filename  in to_be_discarded:
        handle3.write(f"{filename}\n")
~                                                      
