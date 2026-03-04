import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path

from src.data_processors.transcriptomics_processor import TranscriptomicsProcessor
from src.data_processors.proteomics_processor import ProteomicsProcessor
from src.data_processors.genomics_processor import GenomicsProcessor

def process_transcriptomics_data(data_dir: str, cache_dir: str):
    """Process RNA-seq data."""
    print("\n=== Processing Transcriptomics Data ===")
    
    # Initialize processor
    processor = TranscriptomicsProcessor(data_dir=data_dir, cache_dir=cache_dir)
    
    # Example: Load two RNA-seq datasets
    data1 = pd.read_csv(f"{data_dir}/rnaseq_dataset1.csv")
    data2 = pd.read_csv(f"{data_dir}/rnaseq_dataset2.csv")
    
    # Define processing steps
    processing_steps = [
        {
            'name': 'filter_low_expression',
            'type': 'filter_genes',
            'params': {'min_samples': 5}
        },
        {
            'name': 'normalize_data',
            'type': 'normalize',
            'params': {'method': 'log1p', 'scale_factor': 1e6}
        }
    ]
    
    # Process each dataset
    processed1 = processor.process_data(data1, processing_steps)
    processed2 = processor.process_data(data2, processing_steps)
    
    # Calculate quality metrics
    metrics1 = processor.calculate_quality_metrics(processed1)
    print("\nDataset 1 Quality Metrics:")
    print(metrics1)
    
    # Harmonize datasets
    harmonized = processor.harmonize_datasets(
        [processed1, processed2],
        id_types=['ensembl_gene_id', 'symbol'],
        target_type='ensembl_gene_id'
    )
    
    return harmonized

def process_proteomics_data(data_dir: str, cache_dir: str):
    """Process proteomics data."""
    print("\n=== Processing Proteomics Data ===")
    
    # Initialize processor
    processor = ProteomicsProcessor(data_dir=data_dir, cache_dir=cache_dir)
    
    # Example: Load proteomics data
    data = pd.read_csv(f"{data_dir}/proteomics_dataset.csv")
    
    # Define processing steps
    processing_steps = [
        {
            'name': 'filter_proteins',
            'type': 'filter_proteins',
            'params': {'min_samples': 3}
        },
        {
            'name': 'normalize_intensities',
            'type': 'normalize',
            'params': {'method': 'log2'}
        }
    ]
    
    # Process data
    processed = processor.process_data(data, processing_steps)
    
    # Process modifications
    processed = processor.process_modifications(processed, modification_type='phospho')
    
    # Calculate quality metrics
    metrics = processor.calculate_quality_metrics(processed)
    print("\nProteomics Quality Metrics:")
    print(metrics)
    
    # Quantify modifications
    mod_stats = processor.quantify_modifications(processed)
    print("\nModification Statistics:")
    print(mod_stats)
    
    return processed

def process_genomics_data(data_dir: str, cache_dir: str, reference_genome: str):
    """Process genomics data."""
    print("\n=== Processing Genomics Data ===")
    
    # Initialize processor
    processor = GenomicsProcessor(
        data_dir=data_dir,
        cache_dir=cache_dir,
        reference_genome=reference_genome
    )
    
    # Example: Load variant data
    data = pd.read_csv(f"{data_dir}/variants_dataset.csv")
    
    # Define processing steps
    processing_steps = [
        {
            'name': 'filter_variants',
            'type': 'filter_variants',
            'params': {'min_quality': 30, 'min_allele_freq': 0.01}
        },
        {
            'name': 'annotate_variants',
            'type': 'annotate_variants'
        }
    ]
    
    # Process data
    processed = processor.process_data(data, processing_steps)
    
    # Calculate quality metrics
    metrics = processor.calculate_quality_metrics(processed)
    print("\nGenomics Quality Metrics:")
    print(metrics)
    
    return processed

def main():
    # Set up directories
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"
    cache_dir = base_dir / "cache"
    reference_genome = str(data_dir / "reference" / "hg38.fa")
    
    # Create directories if they don't exist
    data_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each data type
    rna_data = process_transcriptomics_data(str(data_dir), str(cache_dir))
    protein_data = process_proteomics_data(str(data_dir), str(cache_dir))
    variant_data = process_genomics_data(str(data_dir), str(cache_dir), reference_genome)
    
    print("\n=== Processing Complete ===")
    print(f"Processed {len(rna_data)} RNA-seq datasets")
    print(f"Processed proteomics data with {len(protein_data)} proteins")
    print(f"Processed genomics data with {len(variant_data)} variants")

if __name__ == "__main__":
    main() 