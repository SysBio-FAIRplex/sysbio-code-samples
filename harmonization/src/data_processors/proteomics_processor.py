import pandas as pd
import numpy as np
import anndata as ad
from typing import Optional, Dict, List, Union
import logging
from pathlib import Path
import scanpy as sc

from .omics_processor import OmicsProcessor

logger = logging.getLogger(__name__)

class ProteomicsProcessor(OmicsProcessor):
    """Processor for proteomics data."""
    
    def __init__(self, data_dir: str, cache_dir: Optional[str] = None):
        """Initialize the proteomics processor."""
        super().__init__(data_dir, cache_dir)
        
    def process_data(self, data: Union[pd.DataFrame, ad.AnnData],
                    processing_steps: List[Dict]) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Process proteomics data according to specified steps.
        
        Args:
            data: Input data
            processing_steps: List of processing steps to apply
            
        Returns:
            Processed data
        """
        processed = data.copy()
        
        for step in processing_steps:
            step_name = step['name']
            step_type = step['type']
            params = step.get('params', {})
            
            logger.info(f"Applying processing step: {step_name}")
            
            try:
                if isinstance(processed, ad.AnnData):
                    if step_type == 'filter_proteins':
                        # Filter proteins by detection frequency
                        min_samples = params.get('min_samples', 3)
                        sc.pp.filter_genes(processed, min_cells=min_samples)
                    elif step_type == 'normalize':
                        method = params.get('method', 'log2')
                        if method == 'log2':
                            processed.X = np.log2(processed.X + 1)
                        elif method == 'zscore':
                            processed.X = (processed.X - np.mean(processed.X, axis=0)) / np.std(processed.X, axis=0)
                    elif step_type == 'batch_correction':
                        sc.pp.combat(processed, key=params.get('batch_key'))
                else:
                    if step_type == 'filter_proteins':
                        min_samples = params.get('min_samples', 3)
                        detection_freq = (processed > 0).sum(axis=0)
                        processed = processed[processed.columns[detection_freq >= min_samples]]
                    elif step_type == 'normalize':
                        method = params.get('method', 'log2')
                        if method == 'log2':
                            processed = np.log2(processed + 1)
                        elif method == 'zscore':
                            processed = (processed - processed.mean()) / processed.std()
                            
            except Exception as e:
                logger.error(f"Error in processing step {step_name}: {str(e)}")
                raise
                
        return processed
        
    def calculate_quality_metrics(self, data: Union[pd.DataFrame, ad.AnnData]) -> Dict:
        """
        Calculate proteomics-specific quality metrics.
        
        Args:
            data: Input data
            
        Returns:
            Dict with quality metrics
        """
        metrics = self.validate_data_quality(data, 'proteomics')
        
        # Add proteomics-specific metrics
        if isinstance(data, ad.AnnData):
            metrics.update({
                'median_proteins_per_sample': np.median(np.sum(data.X > 0, axis=1)),
                'median_intensity': np.median(data.X[data.X > 0]),
                'missing_value_rate': np.mean(data.X == 0),
                'dynamic_range': np.log10(np.percentile(data.X[data.X > 0], 95) / 
                                       np.percentile(data.X[data.X > 0], 5))
            })
        else:
            metrics.update({
                'median_proteins_per_sample': np.median(np.sum(data > 0, axis=1)),
                'median_intensity': np.median(data.values[data.values > 0]),
                'missing_value_rate': np.mean(data == 0),
                'dynamic_range': np.log10(np.percentile(data.values[data.values > 0], 95) / 
                                       np.percentile(data.values[data.values > 0], 5))
            })
            
        return metrics
        
    def process_modifications(self, data: Union[pd.DataFrame, ad.AnnData],
                            modification_type: str) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Process post-translational modifications.
        
        Args:
            data: Input data
            modification_type: Type of modification (e.g., 'phospho', 'ubiq')
            
        Returns:
            Processed data with modification information
        """
        if isinstance(data, ad.AnnData):
            # Add modification information to var annotations
            data.var['modification_type'] = modification_type
            
            # Extract modification sites if available
            if 'protein_site' in data.var:
                data.var['site_position'] = data.var['protein_site'].str.extract('(\d+)').astype(float)
                
        else:
            # Add modification columns
            data['modification_type'] = modification_type
            
            # Extract modification sites if available
            if 'protein_site' in data.columns:
                data['site_position'] = data['protein_site'].str.extract('(\d+)').astype(float)
                
        return data
        
    def quantify_modifications(self, data: Union[pd.DataFrame, ad.AnnData],
                             protein_groups: Optional[pd.DataFrame] = None) -> Dict:
        """
        Quantify post-translational modifications.
        
        Args:
            data: Input data
            protein_groups: Optional protein grouping information
            
        Returns:
            Dict with modification statistics
        """
        stats = {
            'total_sites': 0,
            'modified_proteins': 0,
            'modification_types': {},
            'site_distribution': {}
        }
        
        if isinstance(data, ad.AnnData):
            if 'modification_type' in data.var:
                stats['total_sites'] = data.n_vars
                stats['modified_proteins'] = len(set(data.var.index.str.split('_').str[0]))
                stats['modification_types'] = data.var['modification_type'].value_counts().to_dict()
                
                if 'site_position' in data.var:
                    stats['site_distribution'] = data.var['site_position'].value_counts().to_dict()
        else:
            if 'modification_type' in data.columns:
                stats['total_sites'] = len(data)
                stats['modified_proteins'] = len(set(data.index.str.split('_').str[0]))
                stats['modification_types'] = data['modification_type'].value_counts().to_dict()
                
                if 'site_position' in data.columns:
                    stats['site_distribution'] = data['site_position'].value_counts().to_dict()
                    
        return stats 