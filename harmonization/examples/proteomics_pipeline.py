#!/usr/bin/env python3

import logging
import os
from pathlib import Path

from src.data_processors.proteomics_processor import ProteomicsProcessor
from src.data_processors.amp_pd.proteomics import AMPPDProteomics
from src.data_harmonization.protein_harmonizer import ProteinHarmonizer

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    # Initialize directories
    data_dir = Path("data/proteomics")
    cache_dir = Path("cache/proteomics")
    
    # Create directories if they don't exist
    data_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Initialize processors
    protein_harmonizer = ProteinHarmonizer(cache_dir=cache_dir)
    processor = ProteomicsProcessor(
        data_dir=data_dir,
        cache_dir=cache_dir,
        protein_harmonizer=protein_harmonizer
    )

    # Example processing pipeline for AMP-PD data
    amp_pd_processor = AMPPDProteomics(
        data_dir=data_dir / "amp_pd",
        cache_dir=cache_dir / "amp_pd"
    )

    try:
        # Load data (replace with actual data path)
        raw_data = amp_pd_processor.load_data("example_proteomics.txt")
        
        # Process data
        processed_data = processor.process_data(
            raw_data,
            normalize=True,
            batch_correction=True,
            min_intensity=100,
            missing_value_cutoff=0.3
        )

        # Calculate quality metrics
        qc_metrics = processor.calculate_quality_metrics(processed_data)
        logger.info("Quality metrics: %s", qc_metrics)

        # Process post-translational modifications
        ptm_data = processor.process_modifications(
            processed_data,
            modification_types=["Phosphorylation", "Glycosylation"]
        )

        # Quantify modifications
        mod_stats = processor.quantify_modifications(ptm_data)
        logger.info("Modification statistics: %s", mod_stats)

        # Harmonize protein IDs
        harmonized_data = protein_harmonizer.harmonize_proteins(
            ptm_data,
            source_format="UniProt",
            target_format="RefSeq"
        )

        # Save results
        output_path = data_dir / "processed" / "amp_pd_processed.h5"
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
        logger.error("Error processing proteomics data: %s", str(e))
        raise

if __name__ == "__main__":
    main() 