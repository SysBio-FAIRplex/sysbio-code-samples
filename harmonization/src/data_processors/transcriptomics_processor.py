import pandas as pd
import numpy as np
import anndata as ad
from typing import Optional, Dict, List, Union
import logging
from pathlib import Path
import scanpy as sc

from .omics_processor import OmicsProcessor
from .gene_harmonizer import GeneHarmonizer
from .schema_inferrer import SchemaInferrer

logger = logging.getLogger(__name__)

class TranscriptomicsProcessor(OmicsProcessor):
    """Processor for transcriptomics data."""
    
    def __init__(self, data_dir: str, cache_dir: Optional[str] = None):
        """Initialize the transcriptomics processor."""
        super().__init__(data_dir, cache_dir)
        
        # Initialize harmonizers
        self.gene_harmonizer = GeneHarmonizer(cache_dir=self.cache_dir)
        self.schema_inferrer = SchemaInferrer(model_dir=self.cache_dir)
        
    def load_data(self, file_path: str) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Load transcriptomics data from file.
        
        Args:
            file_path: Path to data file
            
        Returns:
            Loaded data
        """
        file_path = Path(file_path)
        
        try:
            if file_path.suffix == '.h5ad':
                return ad.read_h5ad(file_path)
            elif file_path.suffix in ['.csv', '.txt', '.tsv']:
                return pd.read_csv(file_path, sep='\t' if file_path.suffix == '.tsv' else ',')
            else:
                raise ValueError(f"Unsupported file format: {file_path.suffix}")
                
        except Exception as e:
            logger.error(f"Error loading data from {file_path}: {str(e)}")
            raise
            
    def save_processed_data(self, data: Union[pd.DataFrame, ad.AnnData],
                          output_file: str):
        """
        Save processed data to file.
        
        Args:
            data: Data to save
            output_file: Output file path
        """
        output_path = self.data_dir / output_file
        
        try:
            if isinstance(data, ad.AnnData):
                data.write(output_path)
            else:
                data.to_csv(output_path, sep='\t' if output_path.suffix == '.tsv' else ',')
                
        except Exception as e:
            logger.error(f"Error saving data to {output_path}: {str(e)}")
            raise
            
    def process_data(self, data: Union[pd.DataFrame, ad.AnnData],
                    processing_steps: List[Dict]) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Process transcriptomics data according to specified steps.
        
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
                    if step_type == 'filter_genes':
                        sc.pp.filter_genes(processed, **params)
                    elif step_type == 'normalize':
                        if params.get('method') == 'log1p':
                            sc.pp.normalize_total(processed, target_sum=params.get('scale_factor', 1e6))
                            sc.pp.log1p(processed)
                    elif step_type == 'combat':
                        sc.pp.combat(processed, key=params.get('batch_key'))
                else:
                    if step_type == 'filter_genes':
                        min_cells = params.get('min_cells', 3)
                        gene_counts = (processed > 0).sum()
                        processed = processed[gene_counts[gene_counts >= min_cells].index]
                    elif step_type == 'normalize':
                        if params.get('method') == 'log1p':
                            scale_factor = params.get('scale_factor', 1e6)
                            processed = np.log1p(processed * scale_factor / processed.sum())
                            
            except Exception as e:
                logger.error(f"Error in processing step {step_name}: {str(e)}")
                raise
                
        return processed
        
    def calculate_quality_metrics(self, data: Union[pd.DataFrame, ad.AnnData]) -> Dict:
        """
        Calculate transcriptomics-specific quality metrics.
        
        Args:
            data: Input data
            
        Returns:
            Dict with quality metrics
        """
        metrics = self.validate_data_quality(data, 'transcriptomics')
        
        # Add transcriptomics-specific metrics
        if isinstance(data, ad.AnnData):
            metrics.update({
                'median_genes_per_cell': np.median(np.sum(data.X > 0, axis=1)),
                'median_counts_per_cell': np.median(np.sum(data.X, axis=1)),
                'percent_mito': np.mean(data.var_names.str.startswith('MT-')) if any(data.var_names.str.startswith('MT-')) else 0
            })
        else:
            metrics.update({
                'median_genes_per_sample': np.median(np.sum(data > 0, axis=1)),
                'median_counts_per_sample': np.median(data.sum(axis=1)),
                'percent_mito': np.mean(data.columns.str.startswith('MT-')) if any(data.columns.str.startswith('MT-')) else 0
            })
            
        return metrics
        
    def harmonize_datasets(self, datasets: List[Union[pd.DataFrame, ad.AnnData]],
                         id_types: List[str],
                         target_type: str = 'ensembl_gene_id') -> List[Union[pd.DataFrame, ad.AnnData]]:
        """
        Harmonize multiple datasets.
        
        Args:
            datasets: List of datasets to harmonize
            id_types: List of gene ID types for each dataset
            target_type: Target gene ID type
            
        Returns:
            List of harmonized datasets
        """
        # First infer schemas
        schemas = [self.schema_inferrer.infer_schema(df if isinstance(df, pd.DataFrame)
                                                   else pd.DataFrame(df.X, columns=df.var_names))
                  for df in datasets]
        
        # Get harmonization suggestions
        suggestions = self.schema_inferrer.suggest_harmonization(schemas)
        logger.info(f"Harmonization suggestions: {suggestions}")
        
        # Harmonize gene IDs
        harmonized = self.gene_harmonizer.harmonize_matrices(datasets, id_types, target_type)
        
        # Apply suggested harmonization steps
        for step in suggestions['harmonization_steps']:
            col = step['column']
            action = step['action']
            
            for i, data in enumerate(harmonized):
                if isinstance(data, pd.DataFrame) and col in data.columns:
                    if action == 'convert_to_float':
                        harmonized[i][col] = harmonized[i][col].astype(float)
                    elif action == 'standardize_categories':
                        # Implement category standardization if needed
                        pass
                        
        return harmonized
        
    def get_gene_metadata(self, gene_ids: List[str],
                         id_type: str = 'ensembl_gene_id') -> pd.DataFrame:
        """
        Get metadata for genes.
        
        Args:
            gene_ids: List of gene IDs
            id_type: Type of input IDs
            
        Returns:
            DataFrame with gene metadata
        """
        return self.gene_harmonizer.get_gene_info(gene_ids, id_type) 