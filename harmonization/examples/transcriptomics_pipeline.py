#!/usr/bin/env python3

import logging
import os
from pathlib import Path

from src.data_processors.transcriptomics_processor import TranscriptomicsProcessor
from src.data_processors.amp_ad.transcriptomics import AMPADTranscriptomics
from src.data_harmonization.gene_harmonizer import GeneHarmonizer

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    # Initialize processors
    data_dir = Path("data/transcriptomics")
    cache_dir = Path("cache/transcriptomics")
    
    # Create directories if they don't exist
    data_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Initialize processors
    gene_harmonizer = GeneHarmonizer(cache_dir=cache_dir)
    processor = TranscriptomicsProcessor(
        data_dir=data_dir,
        cache_dir=cache_dir,
        gene_harmonizer=gene_harmonizer
    )

    # Example processing pipeline for AMP-AD data
    amp_ad_processor = AMPADTranscriptomics(
        data_dir=data_dir / "amp_ad",
        cache_dir=cache_dir / "amp_ad"
    )

    try:
        # Load data (replace with actual Synapse IDs)
        raw_data = amp_ad_processor.load_data(synapse_id="syn123456")
        
        # Process data
        processed_data = processor.process_data(
            raw_data,
            normalize=True,
            batch_correction=True,
            min_genes=200,
            min_cells=3,
            percent_mito_cutoff=20
        )

        # Calculate quality metrics
        qc_metrics = processor.calculate_quality_metrics(processed_data)
        logger.info("Quality metrics: %s", qc_metrics)

        # Harmonize gene IDs
        harmonized_data = gene_harmonizer.harmonize_genes(
            processed_data,
            source_format="ENSEMBL",
            target_format="HGNC"
        )

        # Save results
        output_path = data_dir / "processed" / "amp_ad_processed.h5ad"
        processor.save_data(harmonized_data, output_path)
        logger.info("Saved processed data to %s", output_path)

        # Generate quality report
        report_path = data_dir / "processed" / "quality_report.html"
        processor.generate_quality_report(
            harmonized_data,
            qc_metrics,
            report_path
        )
        logger.info("Generated quality report at %s", report_path)

    except Exception as e:
        logger.error("Error processing transcriptomics data: %s", str(e))
        raise

if __name__ == "__main__":
    main() 