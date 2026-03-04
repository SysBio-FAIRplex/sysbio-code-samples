# AD vs PD GWAS Analysis Pipeline

## Overview

This repository contains a comprehensive genome-wide association study (GWAS) pipeline designed to compare Alzheimer's Disease (AD) and Parkinson's Disease (PD) genetic profiles. The analysis performs two main comparisons:

1. **AD Cases vs PD Cases**: Identifies genetic variants that differ between AD and PD patients
2. **AD Controls vs PD Controls**: Identifies genetic variants that differ between healthy controls from AD and PD studies

## Purpose

The primary goal of this analysis is to identify disease-specific genetic associations that may help distinguish between AD and PD at the molecular level. This could provide insights into:

- Disease-specific genetic risk factors
- Potential therapeutic targets
- Biomarkers for differential diagnosis
- Understanding of disease mechanisms

## Data Sources

### AMP-AD (Alzheimer's Disease)
- **Source**: `/path/to/AMP_AD/data/`
- **Cohorts**: ROSMAP, Mayo, MSBB
- **Format**: PLINK pgen/pvar/psam files
- **Population**: European ancestry

### AMP-PD (Parkinson's Disease)
- **Source**: `/path/to/AMP_PD/data/`
- **Format**: PLINK pgen/pvar/psam files
- **Population**: European ancestry

## Pipeline Overview

The analysis consists of 14 main steps:

### 1. Data Preparation
- Clean variant names for consistency
- Standardize chromosome annotations

### 2. Sample Selection
- Extract AD cases and controls from metadata
- Extract PD cases and controls from metadata
- Create sample lists for downstream analysis

### 3. Relatedness Filtering
- Remove related individuals to avoid population stratification
- Process relatedness files from both datasets

### 4. Initial Filtering
- Apply quality control filters:
  - Minor allele frequency (MAF) > 5%
  - Hardy-Weinberg equilibrium (HWE) p > 1e-7
  - SNP-only variants (exclude indels)
  - Missingness < 5%

### 5. Data Merging
- Test initial merges between datasets
- Identify merge conflicts

### 6. Variant Name Resolution
- Resolve merge conflicts by creating positional variant names
- Update variant identifiers for compatibility

### 7. Final Merging
- Successfully merge AD and PD datasets
- Handle palindromic SNPs and other problematic variants

### 8. Phenotype Assignment
- Assign disease labels (AD=2, PD=1)
- Create phenotype files for analysis

### 9. Quality Control
- Apply final quality filters
- Test for missingness patterns
- Remove problematic variants

### 10. Principal Component Analysis
- Calculate principal components for population stratification
- Use thinned datasets for computational efficiency

### 11. Genome-Wide Association Study
- Run logistic regression GWAS
- Include principal components as covariates
- Generate association statistics

### 12. Quality Assessment
- Calculate genomic inflation factor (lambda)
- Assess for population stratification
- Identify significant associations

### 13. Result Filtering
- Remove variants significantly associated in control comparison
- Focus on disease-specific associations

### 14. Allele Frequency Analysis
- Calculate allele frequencies for significant variants
- Prepare results for downstream analysis

## Requirements

### Software Dependencies
- **PLINK 1.9**: For data manipulation and GWAS
- **PLINK 2.0**: For modern file formats and advanced features
- **Python 3.x**: For data processing and analysis
- **Required Python packages**:
  - pandas
  - numpy
  - scipy

### System Requirements
- **Memory**: At least 16GB RAM recommended
- **Storage**: At least 100GB free space
- **CPU**: Multi-core system recommended for parallel processing

### Data Requirements
- Access to AMP-AD and AMP-PD data directories
- Proper permissions for reading genetic data files
- Metadata files for sample annotation

## Installation and Setup

1. **Clone or download this repository**
   ```bash
   git clone [repository-url]
   cd AD_versus_PD_GWAS
   ```

2. **Load required modules** (if using module system)
   ```bash
   module load plink/6-alpha
   module load plink/1.9
   ```

3. **Install Python dependencies**
   ```bash
   pip install pandas numpy scipy
   ```

4. **Update data paths** in the script:
   - Modify `AMP_AD_prefix` and `AMP_PD_prefix` variables
   - Update working directory path

## Usage

### Running the Complete Pipeline

```bash
# Make the script executable
chmod +x AD_versus_PD_GWAS.sh

# Run the complete analysis
./AD_versus_PD_GWAS.sh
```

### Running Individual Steps

The script is modular and can be run step-by-step by commenting out sections. Each step includes progress messages and error checking.

### Monitoring Progress

The script provides detailed progress messages for each step:
- Step completion notifications
- File creation confirmations
- Quality metrics and statistics

## Output Files

### Main Results
- `merged_cases_AD_pheno_GWAS.assoc.logistic`: Case comparison GWAS results
- `merged_controls_AD_pheno_GWAS.assoc.logistic`: Control comparison GWAS results
- `case_versus_case_filtered_results.csv`: Filtered disease-specific results

### Quality Control
- `merged_cases_AD_pheno_GWAS_significant_results.txt`: Significant variants (cases)
- `merged_controls_AD_pheno_GWAS_significant_results.txt`: Significant variants (controls)
- `frequencies_for_case_versus_case.frq`: Allele frequencies

### Intermediate Files
- Various PLINK binary files (.bed, .bim, .fam)
- Principal component files (.eigenvec)
- Sample lists and phenotype files

## Quality Metrics

### Genomic Inflation Factor (Lambda)
- **Cases comparison**: Expected ~1.005 (good control for population stratification)
- **Controls comparison**: Expected ~1.032 (acceptable range)

### Sample Sizes
- **Cases**: 561 AD cases vs 2542 PD cases
- **Controls**: 832 AD controls vs 3031 PD controls

### Quality Filters Applied
- MAF > 5%
- HWE p > 1e-7
- Missingness < 5%
- Relatedness exclusion
- Non-random missingness exclusion

## Interpretation

### Significant Associations
- **P-value threshold**: 5e-8 (genome-wide significance)
- **Effect sizes**: Odds ratios and beta coefficients
- **Direction**: Positive values favor AD, negative values favor PD

### Filtering Strategy
- Variants significantly associated in control comparison are excluded
- Focus on disease-specific genetic differences
- Reduces false positives from population differences

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce thread count or increase system memory
2. **File permission errors**: Check access to data directories
3. **Merge failures**: Check variant naming consistency
4. **Missing dependencies**: Ensure all required software is installed

### Error Recovery
- The script includes intermediate file creation
- Can restart from any step by commenting out completed sections
- Check log files for specific error messages

## Citation

If you use this pipeline in your research, please cite:

```
[Please cite Vitale et 2024 for GenoTools]
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## License

[Completely open source!]

## Contact

For questions or issues, please contact:
- **Email**: [mike@datatecnica.com]
- **Institution**: [DataTecnica]

## Acknowledgments

- AMP-AD consortium for providing AD genetic data
- AMP-PD consortium for providing PD genetic data
- PLINK developers for the genetic analysis software
- Research team and collaborators

---

**Note**: This pipeline is designed for research use and should be validated for clinical applications. Always review results carefully and consider biological plausibility when interpreting findings. 