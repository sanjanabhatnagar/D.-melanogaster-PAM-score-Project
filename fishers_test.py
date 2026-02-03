import scipy.stats as stats

def calculate_enrichment(overlap_bp, total_peak_bp, pam_cluster_bp, total_genome_bp):
    # Creating Contingency Table from Data is in Supplementary Table S3
    table = [
        [overlap_bp, total_peak_bp - overlap_bp],
        [pam_cluster_bp - overlap_bp, total_genome_bp - total_peak_bp - pam_cluster_bp + overlap_bp]
    ]
    odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
    return p_value

# Example - Imp L2 iCLIP dataset
p_val = calculate_enrichment(25, 1555, 38355, 5358482)
print(f"p-value: {p_val:.2e}")
