from abc import ABC, abstractmethod
import pandas as pd
import logging
from pathlib import Path
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class BaseDataProcessor(ABC):
    """Base class for processing data from different sources."""
    
    def __init__(self, data_dir: str):
        """
        Initialize the data processor.
        
        Args:
            data_dir (str): Path to the directory containing the data
        """
        self.data_dir = Path(data_dir)
        self.raw_dir = self.data_dir / 'raw'
        self.processed_dir = self.data_dir / 'processed'
        self._create_directories()
        
    def _create_directories(self):
        """Create necessary directories if they don't exist."""
        os.makedirs(self.raw_dir, exist_ok=True)
        os.makedirs(self.processed_dir, exist_ok=True)
    
    @abstractmethod
    def load_data(self):
        """Load data from source. To be implemented by child classes."""
        pass
    
    @abstractmethod
    def process_data(self):
        """Process the loaded data. To be implemented by child classes."""
        pass
    
    @abstractmethod
    def get_metadata(self):
        """Get metadata about the dataset. To be implemented by child classes."""
        pass
    
    def save_processed_data(self, data: pd.DataFrame, filename: str):
        """
        Save processed data to file.
        
        Args:
            data (pd.DataFrame): Processed data to save
            filename (str): Name of the file to save to
        """
        output_path = self.processed_dir / filename
        data.to_parquet(output_path, index=False)
        logger.info(f"Saved processed data to {output_path}")
    
    def load_processed_data(self, filename: str) -> pd.DataFrame:
        """
        Load previously processed data.
        
        Args:
            filename (str): Name of the file to load
            
        Returns:
            pd.DataFrame: Loaded data
        """
        input_path = self.processed_dir / filename
        if not input_path.exists():
            raise FileNotFoundError(f"No processed data found at {input_path}")
        return pd.read_parquet(input_path)
    
    def validate_data(self, data: pd.DataFrame) -> bool:
        """
        Basic data validation.
        
        Args:
            data (pd.DataFrame): Data to validate
            
        Returns:
            bool: Whether the data is valid
        """
        # Check for empty dataframe
        if data.empty:
            logger.error("Data is empty")
            return False
        
        # Check for missing values
        missing = data.isnull().sum()
        if missing.any():
            logger.warning(f"Missing values found in columns: {missing[missing > 0].to_dict()}")
        
        return True 