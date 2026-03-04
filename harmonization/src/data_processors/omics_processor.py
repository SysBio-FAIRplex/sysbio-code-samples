import pandas as pd
import numpy as np
import anndata as ad
from typing import Optional, Dict, List, Union, Tuple
import logging
from pathlib import Path
import scanpy as sc
from abc import ABC, abstractmethod

from .gene_harmonizer import GeneHarmonizer
from .schema_inferrer import SchemaInferrer

logger = logging.getLogger(__name__)

class OmicsProcessor(ABC):
    """Base class for processing multi-omics data."""
    
    def __init__(self, data_dir: str, cache_dir: Optional[str] = None):
        """
        Initialize the omics processor.
        
        Args:
            data_dir: Directory for data storage
            cache_dir: Directory for caching
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            
        # Initialize harmonizers
        self.gene_harmonizer = GeneHarmonizer(cache_dir=self.cache_dir)
        self.schema_inferrer = SchemaInferrer(model_dir=self.cache_dir)
        
    def load_data(self, file_path: str) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Load omics data from file.
        
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
            elif file_path.suffix == '.parquet':
                return pd.read_parquet(file_path)
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
                if output_path.suffix == '.parquet':
                    data.to_parquet(output_path)
                else:
                    data.to_csv(output_path, sep='\t' if output_path.suffix == '.tsv' else ',')
                
        except Exception as e:
            logger.error(f"Error saving data to {output_path}: {str(e)}")
            raise
            
    @abstractmethod
    def process_data(self, data: Union[pd.DataFrame, ad.AnnData],
                    processing_steps: List[Dict]) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Process omics data according to specified steps.
        Must be implemented by subclasses.
        """
        pass
        
    def infer_and_validate(self, data: Union[pd.DataFrame, ad.AnnData]) -> Tuple[str, Dict]:
        """
        Infer the omics type and validate the data structure.
        
        Args:
            data: Input data
            
        Returns:
            Tuple of (omics_type, validation_results)
        """
        # Convert AnnData to DataFrame for schema inference if needed
        df = data if isinstance(data, pd.DataFrame) else pd.DataFrame(
            data.X, columns=data.var_names, index=data.obs_names
        )
        
        # Infer omics type
        omics_type = self.schema_inferrer.infer_omics_type(df)
        
        # Validate requirements
        validation = self.schema_inferrer.validate_requirements(df, omics_type)
        
        return omics_type, validation
        
    def harmonize_datasets(self, datasets: List[Union[pd.DataFrame, ad.AnnData]],
                         omics_type: Optional[str] = None) -> List[Union[pd.DataFrame, ad.AnnData]]:
        """
        Harmonize multiple datasets.
        
        Args:
            datasets: List of datasets to harmonize
            omics_type: Optional type of omics data
            
        Returns:
            List of harmonized datasets
        """
        # Convert AnnData to DataFrame if needed
        dfs = []
        for data in datasets:
            if isinstance(data, ad.AnnData):
                df = pd.DataFrame(data.X, columns=data.var_names, index=data.obs_names)
                df = pd.concat([df, data.obs], axis=1)
                dfs.append(df)
            else:
                dfs.append(data)
        
        # Infer omics type if not provided
        if omics_type is None:
            omics_types = [self.schema_inferrer.infer_omics_type(df) for df in dfs]
            if len(set(omics_types)) > 1:
                raise ValueError(f"Inconsistent omics types detected: {set(omics_types)}")
            omics_type = omics_types[0]
        
        # Get schemas and harmonization suggestions
        schemas = [self.schema_inferrer.infer_schema(df) for df in dfs]
        suggestions = self.schema_inferrer.suggest_harmonization(schemas, omics_type)
        
        logger.info(f"Harmonization suggestions for {omics_type} data: {suggestions}")
        
        # Apply harmonization steps
        harmonized = []
        for df in dfs:
            processed = df.copy()
            
            # Apply type conversions
            for step in suggestions['harmonization_steps']:
                col = step['column']
                action = step['action']
                
                if col in processed.columns:
                    if action == 'convert_to_float':
                        processed[col] = processed[col].astype(float)
                    elif action == 'standardize_categories':
                        # Implement standardization based on data type
                        if omics_type == 'genomics' and col in ['reference', 'alternate']:
                            processed[col] = processed[col].str.upper()
                
            harmonized.append(processed)
        
        # Convert back to AnnData if original was AnnData
        for i, (data, harm) in enumerate(zip(datasets, harmonized)):
            if isinstance(data, ad.AnnData):
                # Separate expression/abundance data from metadata
                value_cols = [col for col in harm.columns 
                            if harm[col].dtype in [np.float64, np.float32]]
                X = harm[value_cols].values
                
                harmonized[i] = ad.AnnData(
                    X=X,
                    obs=harm.drop(columns=value_cols),
                    var=pd.DataFrame(index=value_cols)
                )
        
        return harmonized
        
    def get_metadata(self, ids: List[str], omics_type: str) -> pd.DataFrame:
        """
        Get metadata for identifiers based on omics type.
        
        Args:
            ids: List of identifiers
            omics_type: Type of omics data
            
        Returns:
            DataFrame with metadata
        """
        if omics_type == 'transcriptomics':
            return self.gene_harmonizer.get_gene_info(ids)
        elif omics_type == 'proteomics':
            # TODO: Implement protein metadata retrieval
            raise NotImplementedError("Protein metadata retrieval not yet implemented")
        elif omics_type == 'genomics':
            # TODO: Implement variant metadata retrieval
            raise NotImplementedError("Variant metadata retrieval not yet implemented")
        else:
            raise ValueError(f"Unknown omics type: {omics_type}")
            
    def validate_data_quality(self, data: Union[pd.DataFrame, ad.AnnData],
                            omics_type: str) -> Dict:
        """
        Perform quality control checks based on omics type.
        
        Args:
            data: Input data
            omics_type: Type of omics data
            
        Returns:
            Dict with quality metrics
        """
        metrics = {
            'sample_count': 0,
            'feature_count': 0,
            'missing_values': 0,
            'warnings': []
        }
        
        if isinstance(data, ad.AnnData):
            metrics['sample_count'] = data.n_obs
            metrics['feature_count'] = data.n_vars
            metrics['missing_values'] = np.isnan(data.X).sum()
        else:
            metrics['sample_count'] = len(data)
            metrics['feature_count'] = len(data.columns)
            metrics['missing_values'] = data.isnull().sum().sum()
        
        # Omics-specific checks
        if omics_type == 'transcriptomics':
            if metrics['feature_count'] < 5000:
                metrics['warnings'].append("Unusually low number of genes detected")
        elif omics_type == 'proteomics':
            if metrics['feature_count'] < 1000:
                metrics['warnings'].append("Unusually low number of proteins detected")
        elif omics_type == 'genomics':
            required_cols = ['chromosome', 'position', 'reference', 'alternate']
            if isinstance(data, pd.DataFrame):
                missing_cols = [col for col in required_cols if col not in data.columns]
                if missing_cols:
                    metrics['warnings'].append(f"Missing required columns: {missing_cols}")
        
        return metrics 