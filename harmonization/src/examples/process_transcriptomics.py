import yaml
from pathlib import Path
import logging
import pandas as pd
import scanpy as sc
import numpy as np

from data_processors.amp_ad.transcriptomics import AMPADTranscriptomics
from harmonization.harmonizer import DataHarmonizer

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Load configuration
    config_path = Path('config/field_mappings.yaml')
    config = load_config(config_path)
    
    # Initialize processors
    amp_ad_processor = AMPADTranscriptomics(
        data_dir='data',
        cache_dir='cache',  # Add cache directory for ML models
        synapse_config='~/.synapseConfig'
    )
    
    # Example Synapse IDs for different brain regions (replace with actual IDs)
    transcriptomics_data = {
        'temporal_cortex': 'syn12345678',
        'hippocampus': 'syn23456789',
        'prefrontal_cortex': 'syn34567890'
    }
    
    try:
        datasets = []
        id_types = []  # Track gene ID types for each dataset
        
        for region, synapse_id in transcriptomics_data.items():
            logger.info(f"Processing {region} data...")
            
            # Load data
            data = amp_ad_processor.load_data_from_synapse(synapse_id)
            
            # Get metadata
            metadata = amp_ad_processor.get_amp_ad_metadata(synapse_id)
            logger.info(f"Dataset metadata: {metadata}")
            
            # Infer schema before processing
            if isinstance(data, pd.DataFrame):
                schema = amp_ad_processor.schema_inferrer.infer_schema(data)
            else:  # AnnData
                schema = amp_ad_processor.schema_inferrer.infer_schema(
                    pd.DataFrame(data.X, columns=data.var_names)
                )
            
            logger.info(f"Inferred schema for {region}:")
            logger.info(f"Number of columns: {schema['metadata']['num_cols']}")
            logger.info(f"Column types: {[info['predicted_type'] for info in schema['columns'].values()]}")
            
            # Determine gene ID type from schema
            gene_cols = [col for col, info in schema['columns'].items() 
                        if info['predicted_type'] == 'gene_id']
            if gene_cols:
                id_type = schema['columns'][gene_cols[0]]['data_type']
                id_types.append(id_type)
                logger.info(f"Detected gene ID type: {id_type}")
            else:
                id_types.append('unknown')
                logger.warning(f"Could not detect gene ID type for {region}")
            
            # Process data with AMP-AD specific steps
            processed_data = amp_ad_processor.process_amp_ad_specific(
                data,
                tissue=region,
                batch_correction=config['rna_seq_parameters']['batch_correction']
            )
            
            # Add region information
            if isinstance(processed_data, pd.DataFrame):
                processed_data['brain_region'] = region
            else:  # AnnData
                processed_data.obs['brain_region'] = region
            
            datasets.append(processed_data)
            
            # Save intermediate results
            amp_ad_processor.save_processed_data(
                processed_data,
                f"{region}_processed.h5ad"
            )
        
        logger.info("Harmonizing datasets...")
        
        # Harmonize datasets using intelligent schema inference
        harmonized_datasets = amp_ad_processor.harmonize_datasets(
            datasets,
            id_types=id_types,
            target_type='ensembl_gene_id'
        )
        
        # Combine harmonized datasets
        if all(isinstance(ds, pd.DataFrame) for ds in harmonized_datasets):
            combined_data = pd.concat(harmonized_datasets, axis=0)
        elif all(isinstance(ds, sc.AnnData) for ds in harmonized_datasets):
            combined_data = sc.concat(harmonized_datasets)
        else:
            raise ValueError("All datasets must be in the same format after harmonization")
        
        # Get metadata for common genes
        if isinstance(combined_data, pd.DataFrame):
            gene_ids = combined_data.columns[combined_data.dtypes == np.float64].tolist()
        else:
            gene_ids = combined_data.var_names.tolist()
            
        gene_metadata = amp_ad_processor.get_gene_metadata(gene_ids)
        logger.info(f"Retrieved metadata for {len(gene_metadata)} genes")
        
        # Save final harmonized data and gene metadata
        amp_ad_processor.save_processed_data(
            combined_data,
            'harmonized_transcriptomics.h5ad'
        )
        gene_metadata.to_csv('data/gene_metadata.tsv', sep='\t', index=False)
        
        logger.info("Processing pipeline completed successfully!")
        
    except Exception as e:
        logger.error(f"Error in transcriptomics processing pipeline: {str(e)}")
        raise

if __name__ == '__main__':
    main() 