import synapseclient
import pandas as pd
from typing import Optional, Dict, List
import logging
from pathlib import Path
import os

from ..base_processor import BaseDataProcessor

logger = logging.getLogger(__name__)

class AMPADProcessor(BaseDataProcessor):
    """Processor for AMP-AD data."""
    
    def __init__(self, data_dir: str, synapse_config: Optional[str] = None):
        """
        Initialize the AMP-AD data processor.
        
        Args:
            data_dir (str): Path to the directory containing the data
            synapse_config (str, optional): Path to synapse configuration file
        """
        super().__init__(data_dir)
        self.syn = synapseclient.Synapse(configPath=synapse_config)
        self.syn.login()
        
    def load_data(self, synapse_id: str) -> pd.DataFrame:
        """
        Load data from Synapse.
        
        Args:
            synapse_id (str): Synapse ID of the dataset
            
        Returns:
            pd.DataFrame: Loaded data
        """
        try:
            entity = self.syn.get(synapse_id)
            if entity.path.endswith('.csv'):
                data = pd.read_csv(entity.path)
            elif entity.path.endswith('.parquet'):
                data = pd.read_parquet(entity.path)
            else:
                raise ValueError(f"Unsupported file format: {entity.path}")
            
            logger.info(f"Loaded data from Synapse ID: {synapse_id}")
            return data
        except Exception as e:
            logger.error(f"Error loading data from Synapse: {str(e)}")
            raise
    
    def process_data(self, data: pd.DataFrame, processing_steps: List[Dict]) -> pd.DataFrame:
        """
        Process the loaded data according to specified steps.
        
        Args:
            data (pd.DataFrame): Data to process
            processing_steps (List[Dict]): List of processing steps to apply
            
        Returns:
            pd.DataFrame: Processed data
        """
        processed_data = data.copy()
        
        for step in processing_steps:
            step_name = step.get('name', 'unknown')
            try:
                if step['type'] == 'filter':
                    processed_data = processed_data.query(step['condition'])
                elif step['type'] == 'transform':
                    processed_data = processed_data.assign(**step['transformations'])
                elif step['type'] == 'aggregate':
                    processed_data = processed_data.groupby(step['by']).agg(step['aggregations']).reset_index()
                
                logger.info(f"Applied processing step: {step_name}")
            except Exception as e:
                logger.error(f"Error in processing step {step_name}: {str(e)}")
                raise
        
        return processed_data
    
    def get_metadata(self, synapse_id: str) -> Dict:
        """
        Get metadata about the dataset.
        
        Args:
            synapse_id (str): Synapse ID of the dataset
            
        Returns:
            Dict: Dataset metadata
        """
        try:
            entity = self.syn.get(synapse_id, downloadFile=False)
            annotations = self.syn.getAnnotations(entity)
            
            metadata = {
                'id': synapse_id,
                'name': entity.name,
                'annotations': annotations,
                'created_on': entity.properties.get('createdOn', ''),
                'modified_on': entity.properties.get('modifiedOn', ''),
                'created_by': entity.properties.get('createdBy', ''),
                'parent_id': entity.properties.get('parentId', '')
            }
            
            return metadata
        except Exception as e:
            logger.error(f"Error getting metadata: {str(e)}")
            raise 