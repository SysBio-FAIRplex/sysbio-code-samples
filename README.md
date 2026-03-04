## sysbio-code-samples

### Description
Code samples for systems biology data analysis, covering genome-wide association studies (GWAS), multi-omics data harmonization, and transcriptomics processing. Built for use with AMP-AD and AMP-PD datasets.

### Usage
Each module is self-contained. See the README in each subdirectory for module-specific usage instructions.

### Data Sources and Attribution
- AMP-AD (Accelerating Medicines Partnership - Alzheimer's Disease)
- AMP-PD (Accelerating Medicines Partnership - Parkinson's Disease)
- License: MIT

### Contributing
Please read [CONTRIBUTING.md](CONTRIBUTING.md) and adhere to our [Code of Conduct](CODE_OF_CONDUCT.md).

### License
MIT License - see [LICENSE](LICENSE) file for details.

### Contact
- Maintainer: Mathew Koretsky (mathew@datatecnica.com)
- Teams: Science and Data Harmonization Working Group, Platform and Tools Working Group

### Acknowledgements
Special thanks to the SysBio FAIRplex Project consortium.

### Project Organization
```
sysbio-code-samples/
├── gwas/
│   ├── README.md
│   ├── AD_versus_PD_GWAS.sh
│   ├── config.yaml
│   └── gwas_analysis_functions.py
├── harmonization/
│   ├── README.md
│   ├── app.py
│   ├── requirements.txt
│   ├── config/
│   │   └── field_mappings.yaml
│   ├── data/
│   │   ├── rnaseq_dataset1.csv
│   │   ├── rnaseq_dataset2.csv
│   │   ├── proteomics_dataset.csv
│   │   └── variants_dataset.csv
│   ├── examples/
│   │   ├── README.md
│   │   ├── transcriptomics_pipeline.py
│   │   ├── proteomics_pipeline.py
│   │   └── multi_omics_example.py
│   └── src/
│       ├── example_usage.py
│       ├── harmonization/
│       │   └── harmonizer.py
│       ├── data_processors/
│       │   ├── base_processor.py
│       │   ├── omics_processor.py
│       │   ├── transcriptomics_processor.py
│       │   ├── proteomics_processor.py
│       │   ├── genomics_processor.py
│       │   ├── gene_harmonizer.py
│       │   ├── schema_inferrer.py
│       │   └── amp_ad/
│       │       ├── processor.py
│       │       └── transcriptomics.py
│       └── examples/
│           └── process_transcriptomics.py
└── transcriptomics/
    └── blood_harmonization.ipynb
```