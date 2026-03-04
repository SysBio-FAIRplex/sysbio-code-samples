#!/bin/bash
# =============================================================================
# AD vs PD GWAS Analysis Pipeline
# =============================================================================
# This script performs a genome-wide association study comparing Alzheimer's 
# Disease (AD) cases vs Parkinson's Disease (PD) cases, and AD controls vs 
# PD controls to identify disease-specific genetic associations.
#
# Author: [Your Name]
# Date: [Date]
# =============================================================================

# =============================================================================
# SETUP AND CONFIGURATION
# =============================================================================

# Set working directory and load required modules
working_dir='/path/to/working/directory/ADvPD'
cd $working_dir
module load plink/6-alpha

# Define data paths
AMP_AD_prefix='/path/to/AMP_AD/data/FILTERED.AMP_AD_EUR_callrate_related_het_geno_haplotype'
AMP_PD_prefix='/path/to/AMP_PD/data/FILTERED.AMP_PD_EUR'

echo "Starting AD vs PD GWAS analysis..."
echo "Working directory: $working_dir"

# =============================================================================
# STEP 1: DATA PREPARATION - Clean variant names
# =============================================================================
echo "Step 1: Cleaning variant names..."

# Clean up variant names for AD data
awk '{print $3 "\t" "chr" $3}' $AMP_AD_prefix.pvar > AMP_AD_variant_names_to_update.txt
plink2 --pfile $AMP_AD_prefix --update-name AMP_AD_variant_names_to_update.txt --make-pgen --out AMP_AD_raw

echo "Variant names cleaned for AD data"

# =============================================================================
# STEP 2: SAMPLE SELECTION - Extract cases and controls
# =============================================================================
echo "Step 2: Extracting cases and controls..."

# Create Python script for sample selection
cat > extract_samples.py << 'EOF'
import pandas as pd

# Extract PD samples
print("Processing PD samples...")
amp_pd_meta_df = pd.read_csv("/path/to/AMP_PD/metadata/AMPPD_v3_COV_wPHENOS_wGENOTOOLS.csv")
amp_pd_control_subset_df = amp_pd_meta_df[amp_pd_meta_df['LATEST_DX'] == 'No PD Nor Other Neurological Disorder']
amp_pd_case_subset_df = amp_pd_meta_df[amp_pd_meta_df['LATEST_DX'].isin(['Parkinson\'s Disease', 'Idiopathic PD'])]

# Save PD sample lists
amp_pd_control_subset_df[["ID","ID"]].to_csv("amp_pd_control_list", sep="\t", index=False)
amp_pd_case_subset_df[["ID","ID"]].to_csv("amp_pd_case_list", sep="\t", index=False)

# Extract AD samples
print("Processing AD samples...")
amp_ad_meta_df = pd.read_csv("/path/to/AMP_AD/metadata/threeCohorts.csv")
amp_ad_meta_df['ID'] = amp_ad_meta_df['wgs_id'].astype(str).str.split('.').str[0]
amp_ad_control_subset_df = amp_ad_meta_df[amp_ad_meta_df['clinAD'] == 0]
amp_ad_case_subset_df = amp_ad_meta_df[amp_ad_meta_df['clinAD'] == 1]

# Save AD sample lists
amp_ad_control_subset_df[["ID","ID"]].to_csv("amp_ad_control_list", sep="\t", index=False)
amp_ad_case_subset_df[["ID","ID"]].to_csv("amp_ad_case_list", sep="\t", index=False)

print("Sample lists created successfully")
EOF

# Run the Python script
python extract_samples.py

# =============================================================================
# STEP 3: RELATEDNESS FILTERING - Remove related individuals
# =============================================================================
echo "Step 3: Processing relatedness files..."

# Check relatedness file structure
echo "Checking PD relatedness file structure..."
echo "Column 1 unique values:" && awk -F',' '{print $1}' /path/to/AMP_PD/related_files/FILTERED.AMP_PD_EUR.related | sort | uniq | wc -l && echo "Column 3 unique values:" && awk -F',' '{print $3}' /path/to/AMP_PD/related_files/FILTERED.AMP_PD_EUR.related | sort | uniq | wc -l

echo "Checking AD relatedness file structure..."
echo "Column 1 unique values:" && awk -F',' '{print $1}' /path/to/AMP_AD/genotools/FILTERED.AMP_AD_EUR.related | sort | uniq | wc -l && echo "Column 3 unique values:" && awk -F',' '{print $3}' /path/to/AMP_AD/genotools/FILTERED.AMP_AD_EUR.related | sort | uniq | wc -l

# Create exclusion lists for related individuals
awk -F',' 'NR==1{print $1"\t"$2; next} {print $1"\t"$2}' /path/to/AMP_PD/related_files/FILTERED.AMP_PD_EUR.related > amp_pd_related_to_exclude.tsv
awk -F',' 'NR==1{print $1"\t"$2; next} {print $1"\t"$2}' /path/to/AMP_AD/genotools/FILTERED.AMP_AD_EUR.related > amp_ad_related_to_exclude.tsv

# =============================================================================
# STEP 4: INITIAL FILTERING - Create case/control subsets with quality filters
# =============================================================================
echo "Step 4: Creating filtered case/control subsets..."

# Filter AD samples
echo "Filtering AD cases..."
plink2 --pfile AMP_AD_raw --keep amp_ad_case_list --remove amp_ad_related_to_exclude.tsv \
    --maf 0.05 --hwe 0.0000001 --snps-only 'just-acgt' --geno 0.05 \
    --make-bed --out amp_ad_cases

echo "Filtering AD controls..."
plink2 --pfile AMP_AD_raw --keep amp_ad_control_list --remove amp_ad_related_to_exclude.tsv \
    --maf 0.05 --hwe 0.0000001 --snps-only 'just-acgt' --geno 0.05 \
    --make-bed --out amp_ad_controls

# Filter PD samples
echo "Filtering PD cases..."
plink2 --pfile $AMP_PD_prefix --keep amp_pd_case_list --remove amp_pd_related_to_exclude.tsv \
    --maf 0.05 --hwe 0.0000001 --snps-only 'just-acgt' --geno 0.05 \
    --make-bed --out amp_pd_cases

echo "Filtering PD controls..."
plink2 --pfile $AMP_PD_prefix --keep amp_pd_control_list --remove amp_pd_related_to_exclude.tsv \
    --maf 0.05 --hwe 0.0000001 --snps-only 'just-acgt' --geno 0.05 \
    --make-bed --out amp_pd_controls

# =============================================================================
# STEP 5: DATA MERGING - Test initial merges
# =============================================================================
echo "Step 5: Testing initial data merges..."

module load plink/1.9

# Test merging controls
echo "Testing control merge..."
plink --bfile amp_pd_controls --bmerge amp_ad_controls.bed amp_ad_controls.bim amp_ad_controls.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_controls

# Test merging cases
echo "Testing case merge..."
plink --bfile amp_pd_cases --bmerge amp_ad_cases.bed amp_ad_cases.bim amp_ad_cases.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_cases

# =============================================================================
# STEP 6: VARIANT NAME RESOLUTION - Handle merge conflicts
# =============================================================================
echo "Step 6: Resolving variant name conflicts..."

# Create positional variant names to resolve merge conflicts
awk '{print $2"\t""chr"$1":"$4}' amp_ad_cases.bim > amp_ad_cases.positional_varname_update.tsv
awk '{print $2"\t""chr"$1":"$4}' amp_ad_controls.bim > amp_ad_controls.positional_varname_update.tsv
awk '{print $2"\t""chr"$1":"$4}' amp_pd_cases.bim > amp_pd_cases.positional_varname_update.tsv
awk '{print $2"\t""chr"$1":"$4}' amp_pd_controls.bim > amp_pd_controls.positional_varname_update.tsv

# Update variant names using PLINK2
module load plink/6-alpha

echo "Updating variant names for cases..."
plink2 --bfile amp_ad_cases --update-name amp_ad_cases.positional_varname_update.tsv --make-bed --out amp_ad_cases_positional_varnames
plink2 --bfile amp_pd_cases --update-name amp_pd_cases.positional_varname_update.tsv --make-bed --out amp_pd_cases_positional_varnames

echo "Updating variant names for controls..."
plink2 --bfile amp_ad_controls --update-name amp_ad_controls.positional_varname_update.tsv --make-bed --out amp_ad_controls_positional_varnames
plink2 --bfile amp_pd_controls --update-name amp_pd_controls.positional_varname_update.tsv --make-bed --out amp_pd_controls_positional_varnames

# =============================================================================
# STEP 7: FINAL MERGING - Merge with updated variant names
# =============================================================================
echo "Step 7: Final data merging..."

module load plink/1.9

# Test merges with updated variant names
echo "Testing control merge with updated names..."
plink --bfile amp_pd_controls_positional_varnames --bmerge amp_ad_controls_positional_varnames.bed amp_ad_controls_positional_varnames.bim amp_ad_controls_positional_varnames.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_controls_positional_varnames_test

echo "Testing case merge with updated names..."
plink --bfile amp_pd_cases_positional_varnames --bmerge amp_ad_cases_positional_varnames.bed amp_ad_cases_positional_varnames.bim amp_ad_cases_positional_varnames.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_cases_positional_varnames_test

# Remove problematic variants from base datasets
echo "Removing problematic variants..."
plink --bfile amp_pd_controls_positional_varnames --exclude merged_controls_positional_varnames_test-merge.missnp --make-bed --out amp_pd_controls_positional_varnames_to_merge
plink --bfile amp_pd_cases_positional_varnames --exclude merged_cases_positional_varnames_test-merge.missnp --make-bed --out amp_pd_cases_positional_varnames_to_merge

# Final merge attempts
echo "Final control merge..."
plink --bfile amp_pd_controls_positional_varnames_to_merge --bmerge amp_ad_controls_positional_varnames.bed amp_ad_controls_positional_varnames.bim amp_ad_controls_positional_varnames.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_controls_positional_varnames_test_next

echo "Final case merge..."
plink --bfile amp_pd_cases_positional_varnames_to_merge --bmerge amp_ad_cases_positional_varnames.bed amp_ad_cases_positional_varnames.bim amp_ad_cases_positional_varnames.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_cases_positional_varnames_test_next

# Handle remaining problematic variants (palindromic SNPs)
echo "Handling palindromic SNPs..."
plink --bfile amp_pd_controls_positional_varnames_to_merge --exclude merged_controls_positional_varnames_test-merge.missnp --make-bed --out amp_pd_controls_positional_varnames_to_merge_pals
plink --bfile amp_pd_cases_positional_varnames_to_merge --exclude merged_cases_positional_varnames_test-merge.missnp --make-bed --out amp_pd_cases_positional_varnames_to_merge_pals
plink --bfile amp_ad_controls_positional_varnames --exclude merged_controls_positional_varnames_test-merge.missnp --make-bed --out amp_ad_controls_positional_varnames_to_merge_pals
plink --bfile amp_ad_cases_positional_varnames --exclude merged_cases_positional_varnames_test-merge.missnp --make-bed --out amp_ad_cases_positional_varnames_to_merge_pals

# Final successful merges
echo "Final successful merges..."
plink --bfile amp_pd_controls_positional_varnames_to_merge_pals --bmerge amp_ad_controls_positional_varnames_to_merge_pals.bed amp_ad_controls_positional_varnames_to_merge_pals.bim amp_ad_controls_positional_varnames_to_merge_pals.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_controls_positional_varnames_test_next

plink --bfile amp_pd_cases_positional_varnames_to_merge_pals --bmerge amp_ad_cases_positional_varnames_to_merge_pals.bed amp_ad_cases_positional_varnames_to_merge_pals.bim amp_ad_cases_positional_varnames_to_merge_pals.fam \
    --geno 0.05 --hwe 0.00001 --make-bed --out merged_cases_positional_varnames_test_next

# =============================================================================
# STEP 8: PHENOTYPE ASSIGNMENT - Add disease labels
# =============================================================================
echo "Step 8: Assigning phenotypes..."

# Create phenotype files (AD=2, PD=1)
awk '{print $1"\t"$2"\t""2"}' amp_ad_controls_positional_varnames_to_merge_pals.fam > temp_ad_control.pheno
awk '{print $1"\t"$2"\t""2"}' amp_ad_cases_positional_varnames_to_merge_pals.fam > temp_ad_case.pheno
awk '{print $1"\t"$2"\t""1"}' amp_pd_controls_positional_varnames_to_merge_pals.fam > temp_pd_control.pheno
awk '{print $1"\t"$2"\t""1"}' amp_pd_cases_positional_varnames_to_merge_pals.fam > temp_pd_case.pheno

# Combine phenotype files
cat temp_ad_case.pheno temp_pd_case.pheno > all_cases.pheno
cat temp_ad_control.pheno temp_pd_control.pheno > all_controls.pheno

# Add phenotypes to merged datasets
plink --bfile merged_cases_positional_varnames_test_next --pheno all_cases.pheno --make-bed --out merged_cases_AD_pheno
plink --bfile merged_controls_positional_varnames_test_next --pheno all_controls.pheno --make-bed --out merged_controls_AD_pheno

# =============================================================================
# STEP 9: QUALITY CONTROL - Final filtering and quality checks
# =============================================================================
echo "Step 9: Final quality control..."

# Apply final quality filters
echo "Applying MAF, HWE, and missingness filters..."
plink --bfile merged_cases_AD_pheno --maf 0.05 --hwe 0.00001 --geno 0.05 --make-bed --out merged_cases_AD_pheno_maf_geno_hwe --allow-no-sex
plink --bfile merged_controls_AD_pheno --maf 0.05 --hwe 0.00001 --geno 0.05 --make-bed --out merged_controls_AD_pheno_maf_geno_hwe --allow-no-sex

# Test for missingness by haplotype and non-random missingness
echo "Testing missingness patterns..."
plink --bfile merged_cases_AD_pheno_maf_geno_hwe --test-mishap --out merged_cases_AD_pheno_maf_geno_hwe_misshap_result --allow-no-sex
plink --bfile merged_cases_AD_pheno_maf_geno_hwe --test-missing --out merged_cases_AD_pheno_maf_geno_hwe_testmissing_result --allow-no-sex

plink --bfile merged_controls_AD_pheno_maf_geno_hwe --test-mishap --out merged_controls_AD_pheno_maf_geno_hwe_misshap_result --allow-no-sex
plink --bfile merged_controls_AD_pheno_maf_geno_hwe --test-missing --out merged_controls_AD_pheno_maf_geno_hwe_testmissing_result --allow-no-sex

# Create exclusion list for problematic variants
echo "Creating exclusion list for problematic variants..."
cat merged_cases_AD_pheno_maf_geno_hwe_testmissing_result.missing merged_controls_AD_pheno_maf_geno_hwe_testmissing_result.missing | awk '$5 < 0.01 {print $2"\t"$5}' > test_missing_to_exclude.txt
cat merged_cases_AD_pheno_maf_geno_hwe_misshap_result.missing.hap merged_controls_AD_pheno_maf_geno_hwe_misshap_result.missing.hap | awk '$8 < 0.01 {print $1"\t"$8}' > mishap_to_exclude.txt

cat test_missing_to_exclude.txt mishap_to_exclude.txt > mishap_and_testmiss_exclusion_list.txt

# Apply exclusions to create final analysis datasets
echo "Creating final analysis-ready datasets..."
plink --bfile merged_controls_AD_pheno_maf_geno_hwe --exclude mishap_and_testmiss_exclusion_list.txt --make-bed --out merged_controls_AD_pheno_analysis_ready
plink --bfile merged_cases_AD_pheno_maf_geno_hwe --exclude mishap_and_testmiss_exclusion_list.txt --make-bed --out merged_cases_AD_pheno_analysis_ready

# =============================================================================
# STEP 10: PRINCIPAL COMPONENT ANALYSIS - Population stratification
# =============================================================================
echo "Step 10: Principal component analysis..."

# Thin data for PCA (reduce computational burden)
echo "Thinning data for PCA..."
plink --bfile merged_controls_AD_pheno_analysis_ready --thin 0.001 --make-bed --out merged_controls_AD_pheno_analysis_ready_thinned --allow-no-sex
plink --bfile merged_cases_AD_pheno_analysis_ready --thin 0.001 --make-bed --out merged_cases_AD_pheno_analysis_ready_thinned --allow-no-sex

# Calculate principal components
echo "Calculating principal components..."
plink --bfile merged_controls_AD_pheno_analysis_ready_thinned --pca 10 header tabs --out merged_controls_AD_pheno_analysis_ready_PCs --allow-no-sex
plink --bfile merged_cases_AD_pheno_analysis_ready_thinned --pca 10 header tabs --out merged_cases_AD_pheno_analysis_ready_PCs --allow-no-sex

# =============================================================================
# STEP 11: GENOME-WIDE ASSOCIATION STUDY - Main analysis
# =============================================================================
echo "Step 11: Running genome-wide association studies..."

# Run GWAS for controls comparison (AD controls vs PD controls)
echo "Running GWAS for controls comparison..."
plink --bfile merged_controls_AD_pheno_analysis_ready --covar merged_controls_AD_pheno_analysis_ready_PCs.eigenvec \
    --out merged_controls_AD_pheno_GWAS --allow-no-sex --logistic beta --threads 10

# Run GWAS for cases comparison (AD cases vs PD cases)
echo "Running GWAS for cases comparison..."
plink --bfile merged_cases_AD_pheno_analysis_ready --covar merged_cases_AD_pheno_analysis_ready_PCs.eigenvec \
    --out merged_cases_AD_pheno_GWAS --allow-no-sex --logistic beta --threads 10

# =============================================================================
# STEP 12: QUALITY ASSESSMENT - Lambda calculation and QQ plots
# =============================================================================
echo "Step 12: Quality assessment and lambda calculation..."

# Create Python script for quality assessment
cat > quality_assessment.py << 'EOF'
import pandas as pd
import numpy as np
from scipy import stats

print("Analyzing GWAS results...")

# Analyze cases comparison results
print("\n=== Cases Comparison (AD cases vs PD cases) ===")
gwas_results = pd.read_csv('merged_cases_AD_pheno_GWAS.assoc.logistic', delim_whitespace=True, engine='c')
gwas_results = gwas_results[gwas_results['TEST'] == "ADD"]

# Convert P-values to chi-square statistics
gwas_results['CHISQ'] = stats.chi2.ppf(1 - gwas_results['P'], 1)

# Calculate lambda
chi2 = gwas_results['CHISQ'].dropna()
lambda_gc = np.median(chi2) / 0.455  # 0.455 is the median of chi2 distribution with 1 df

# Calculate scaled lambda for 1000 cases and 1000 controls
n_cases = 561
n_controls = 2542
lambda_scaled = 1 + (lambda_gc - 1) * ((1/n_cases + 1/n_controls) / (1/1000 + 1/1000))

print(f"Lambda: {lambda_gc:.4f}")
print(f"Scaled Lambda: {lambda_scaled:.4f}")

# Extract significant results (P < 5e-8)
significant_results = gwas_results[gwas_results['P'] < 0.00000005]
significant_results.to_csv('merged_cases_AD_pheno_GWAS_significant_results.txt', sep='\t', index=False)
print(f"Number of significant variants: {len(significant_results)}")

# Analyze controls comparison results
print("\n=== Controls Comparison (AD controls vs PD controls) ===")
gwas_results = pd.read_csv('merged_controls_AD_pheno_GWAS.assoc.logistic', delim_whitespace=True, engine='c')
gwas_results = gwas_results[gwas_results['TEST'] == "ADD"]

# Convert P-values to chi-square statistics
gwas_results['CHISQ'] = stats.chi2.ppf(1 - gwas_results['P'], 1)

# Calculate lambda
chi2 = gwas_results['CHISQ'].dropna()
lambda_gc = np.median(chi2) / 0.455

# Calculate scaled lambda for 1000 cases and 1000 controls
n_cases = 832
n_controls = 3031
lambda_scaled = 1 + (lambda_gc - 1) * ((1/n_cases + 1/n_controls) / (1/1000 + 1/1000))

print(f"Lambda: {lambda_gc:.4f}")
print(f"Scaled Lambda: {lambda_scaled:.4f}")

# Extract significant results (P < 5e-8)
significant_results = gwas_results[gwas_results['P'] < 0.00000005]
significant_results.to_csv('merged_controls_AD_pheno_GWAS_significant_results.txt', sep='\t', index=False)
print(f"Number of significant variants: {len(significant_results)}")

print("\nQuality assessment complete!")
EOF

# Run quality assessment
python quality_assessment.py

# =============================================================================
# STEP 13: RESULT FILTERING - Remove control-associated variants from case analysis
# =============================================================================
echo "Step 13: Filtering results to remove control-associated variants..."

# Create Python script for result filtering
cat > filter_results.py << 'EOF'
import pandas as pd

print("Filtering case results to remove control-associated variants...")

# Read GWAS results
cases_df = pd.read_csv('merged_cases_AD_pheno_GWAS.assoc.logistic', delim_whitespace=True, engine='c')
controls_df = pd.read_csv('merged_controls_AD_pheno_GWAS.assoc.logistic', delim_whitespace=True, engine='c')

# Filter for additive model results only
filtered_cases = cases_df[cases_df['TEST'].isin(["ADD"])]
filtered_controls = controls_df[controls_df['TEST'].isin(["ADD"])]

# Get SNPs from controls where P > 1E-5 (not significantly associated in controls)
insignificant_control_snps = filtered_controls[filtered_controls['P'] > 1E-5]['SNP']

# Filter cases for those SNPs (remove variants significantly associated in controls)
final_df = filtered_cases[filtered_cases['SNP'].isin(insignificant_control_snps)]

# Export filtered results
final_df.to_csv('case_versus_case_filtered_results.csv', index=False)
print(f"Filtered results saved. Remaining variants: {len(final_df)}")
EOF

# Run result filtering
python filter_results.py

# =============================================================================
# STEP 14: ALLELE FREQUENCY ANALYSIS - Calculate frequencies for significant variants
# =============================================================================
echo "Step 14: Calculating allele frequencies for significant variants..."

# Calculate allele frequencies for the case comparison dataset
plink --bfile merged_cases_AD_pheno_analysis_ready --freq --out frequencies_for_case_versus_case

echo "Allele frequency calculation complete!"

# =============================================================================
# COMPLETION
# =============================================================================
echo "=========================================="
echo "AD vs PD GWAS Analysis Complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "- merged_cases_AD_pheno_GWAS.assoc.logistic: Case comparison GWAS results"
echo "- merged_controls_AD_pheno_GWAS.assoc.logistic: Control comparison GWAS results"
echo "- case_versus_case_filtered_results.csv: Filtered case results"
echo "- frequencies_for_case_versus_case.frq: Allele frequencies"
echo "- Quality assessment results printed above"
echo ""
echo "Analysis completed successfully!" 