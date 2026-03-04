import pandas as pd
import numpy as np
from typing import Dict, List, Optional
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class DataHarmonizer:
    """Class for harmonizing data across different sources."""
    
    def __init__(self, mapping_config: Dict):
        """
        Initialize the harmonizer with a mapping configuration.
        
        Args:
            mapping_config (Dict): Configuration for mapping fields across datasets
        """
        self.mapping_config = mapping_config
        
    def harmonize_columns(self, data: pd.DataFrame, source: str) -> pd.DataFrame:
        """
        Harmonize column names based on mapping configuration.
        
        Args:
            data (pd.DataFrame): Data to harmonize
            source (str): Source of the data (e.g., 'AMP-AD', 'AMP-PD', 'CMDGA')
            
        Returns:
            pd.DataFrame: Data with harmonized column names
        """
        if source not in self.mapping_config:
            raise ValueError(f"No mapping configuration found for source: {source}")
            
        mapping = self.mapping_config[source]
        return data.rename(columns=mapping)
    
    def standardize_values(self, data: pd.DataFrame, column: str, 
                         standard_values: Dict[str, str]) -> pd.DataFrame:
        """
        Standardize values in a column according to a mapping.
        
        Args:
            data (pd.DataFrame): Data to standardize
            column (str): Column to standardize
            standard_values (Dict[str, str]): Mapping of original values to standard values
            
        Returns:
            pd.DataFrame: Data with standardized values
        """
        if column not in data.columns:
            logger.warning(f"Column {column} not found in data")
            return data
            
        standardized = data.copy()
        standardized[column] = standardized[column].map(standard_values).fillna(standardized[column])
        return standardized
    
    def harmonize_units(self, data: pd.DataFrame, column: str, 
                       current_unit: str, target_unit: str, 
                       conversion_factor: float) -> pd.DataFrame:
        """
        Convert values from one unit to another.
        
        Args:
            data (pd.DataFrame): Data to convert
            column (str): Column to convert
            current_unit (str): Current unit of measurement
            target_unit (str): Target unit of measurement
            conversion_factor (float): Factor to multiply values by
            
        Returns:
            pd.DataFrame: Data with converted values
        """
        if column not in data.columns:
            logger.warning(f"Column {column} not found in data")
            return data
            
        converted = data.copy()
        converted[column] = converted[column] * conversion_factor
        converted[f"{column}_unit"] = target_unit
        return converted
    
    def merge_datasets(self, datasets: List[pd.DataFrame], 
                      on: List[str], how: str = 'outer') -> pd.DataFrame:
        """
        Merge multiple harmonized datasets.
        
        Args:
            datasets (List[pd.DataFrame]): List of datasets to merge
            on (List[str]): Columns to merge on
            how (str): Type of merge to perform
            
        Returns:
            pd.DataFrame: Merged dataset
        """
        if not datasets:
            raise ValueError("No datasets provided for merging")
            
        result = datasets[0]
        for i, dataset in enumerate(datasets[1:], 1):
            try:
                result = result.merge(dataset, on=on, how=how, 
                                    suffixes=(f'_1', f'_{i+1}'))
            except Exception as e:
                logger.error(f"Error merging dataset {i+1}: {str(e)}")
                raise
                
        return result
    
    def validate_harmonization(self, original_data: pd.DataFrame, 
                             harmonized_data: pd.DataFrame, 
                             expected_columns: List[str]) -> bool:
        """
        Validate the harmonization process.
        
        Args:
            original_data (pd.DataFrame): Original data
            harmonized_data (pd.DataFrame): Harmonized data
            expected_columns (List[str]): Expected columns after harmonization
            
        Returns:
            bool: Whether the harmonization is valid
        """
        # Check if all expected columns are present
        missing_columns = set(expected_columns) - set(harmonized_data.columns)
        if missing_columns:
            logger.error(f"Missing expected columns: {missing_columns}")
            return False
            
        # Check if row count is preserved
        if len(original_data) != len(harmonized_data):
            logger.error("Row count changed during harmonization")
            return False
            
        # Check for unexpected null values
        new_nulls = harmonized_data[expected_columns].isnull().sum()
        if new_nulls.any():
            logger.warning(f"New null values introduced in columns: {new_nulls[new_nulls > 0].to_dict()}")
            
        return True 