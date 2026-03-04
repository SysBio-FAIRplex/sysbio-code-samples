import yaml
from pathlib import Path
import logging

from data_processors.amp_ad.processor import AMPADProcessor
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
    
    # Initialize processors and harmonizer
    amp_ad_processor = AMPADProcessor(
        data_dir='data',
        synapse_config='~/.synapseConfig'  # Update with your Synapse config path
    )
    
    harmonizer = DataHarmonizer(mapping_config=config['common_fields'])
    
    # Example Synapse IDs (replace with actual IDs)
    clinical_data_id = "syn12345678"
    
    try:
        # Load data
        logger.info("Loading AMP-AD clinical data...")
        clinical_data = amp_ad_processor.load_data(clinical_data_id)
        
        # Get metadata
        metadata = amp_ad_processor.get_metadata(clinical_data_id)
        logger.info(f"Dataset metadata: {metadata}")
        
        # Process data
        processing_steps = [
            {
                'name': 'filter_complete_cases',
                'type': 'filter',
                'condition': 'age.notna() and diagnosis.notna()'
            },
            {
                'name': 'add_source',
                'type': 'transform',
                'transformations': {
                    'study_source': '"AMP-AD"'
                }
            }
        ]
        
        processed_data = amp_ad_processor.process_data(clinical_data, processing_steps)
        
        # Harmonize data
        logger.info("Harmonizing data...")
        harmonized_data = harmonizer.harmonize_columns(processed_data, source='AMP-AD')
        
        # Standardize sex values
        harmonized_data = harmonizer.standardize_values(
            harmonized_data,
            column='sex',
            standard_values=config['sex_standardization']
        )
        
        # Validate harmonization
        is_valid = harmonizer.validate_harmonization(
            processed_data,
            harmonized_data,
            expected_columns=config['expected_columns']
        )
        
        if is_valid:
            logger.info("Data harmonization successful!")
            # Save harmonized data
            amp_ad_processor.save_processed_data(
                harmonized_data,
                'harmonized_clinical_data.parquet'
            )
        else:
            logger.error("Data harmonization failed validation!")
            
    except Exception as e:
        logger.error(f"Error in data processing pipeline: {str(e)}")
        raise

if __name__ == '__main__':
    main() 