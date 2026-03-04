# Example Processing Pipelines

This directory contains example processing pipelines for different types of omics data using our data harmonization framework.

## Available Examples

1. `transcriptomics_pipeline.py`: Example pipeline for processing RNA-seq data from AMP-AD
   - Demonstrates data loading, quality control, normalization, and batch correction
   - Shows gene ID harmonization between different formats
   - Includes quality metrics calculation and reporting

2. `proteomics_pipeline.py`: Example pipeline for processing proteomics data from AMP-PD
   - Shows protein data processing with intensity normalization
   - Handles post-translational modifications
   - Includes protein ID harmonization
   - Demonstrates quality control and reporting

3. `genomics_pipeline.py`: Example pipeline for processing genomics data from CMDGA
   - Demonstrates variant processing and filtering
   - Shows coordinate conversion between genome builds
   - Includes variant annotation and merging
   - Calculates various genomics-specific quality metrics

## Prerequisites

Before running these examples, make sure you have:

1. Installed all required dependencies (see main README.md)
2. Set up appropriate data directories
3. Configured access to required data sources (e.g., Synapse credentials for AMP-AD data)
4. Downloaded necessary reference files (e.g., genome reference for variant annotation)

## Directory Structure

The examples expect the following directory structure:

```
project_root/
├── data/
│   ├── transcriptomics/
│   ├── proteomics/
│   ├── genomics/
│   └── reference/
│       └── genome/
├── cache/
│   ├── transcriptomics/
│   ├── proteomics/
│   └── genomics/
└── examples/
    ├── transcriptomics_pipeline.py
    ├── proteomics_pipeline.py
    └── genomics_pipeline.py
```

## Running the Examples

1. First, create the necessary directories:
```bash
mkdir -p data/{transcriptomics,proteomics,genomics,reference/genome}
mkdir -p cache/{transcriptomics,proteomics,genomics}
```

2. Download or prepare your input data files and place them in the appropriate directories.

3. Run an example pipeline:
```bash
python examples/transcriptomics_pipeline.py
```

## Customizing the Examples

Each example pipeline is designed to be easily customizable:

- Modify processing parameters in the `process_data()` calls
- Adjust quality control thresholds
- Change input/output file paths
- Add or remove processing steps
- Modify harmonization settings

## Output

Each pipeline generates:

1. Processed data files in standard formats (e.g., h5ad, h5, vcf.gz)
2. Quality control reports in HTML format
3. Log files with processing information and statistics

## Notes

- These are example pipelines and may need to be modified based on your specific data and requirements
- Make sure to replace placeholder values (e.g., Synapse IDs, file paths) with your actual values
- Consider memory requirements when processing large datasets
- Some processing steps may take significant time to complete 