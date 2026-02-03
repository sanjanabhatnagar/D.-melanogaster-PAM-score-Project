import scipy.stats as stats

def calculate_enrichment(overlap_bp, total_peak_bp, pam_group_bp, total_genome_bp):
    # Creating Contingency Table from Data is in Supplementary Table S3
    table = [
        [overlap_bp, total_peak_bp - overlap_bp], # Is a peak-  overlap_bp = Inside a PAM group and (total_peak_bp - overlap_bp) = Outside the PAM group
        [pam_group_bp - overlap_bp, total_genome_bp - total_peak_bp - pam_group_bp + overlap_bp] # Is not a peak -  pam_cluster_bp - overlap_bp = PAM group no binding and 
        # (total_genome_bp - total_peak_bp - pam_group_bp + overlap_bp) = the empty background - no peaks
    ]
    odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
    return p_value

# Example - Imp L2 iCLIP dataset
p_val = calculate_enrichment(25, 1555, 38355, 5358482)
print(f"p-value: {p_val:.2e}")
