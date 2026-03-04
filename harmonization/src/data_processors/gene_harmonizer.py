from typing import Dict, List, Set, Optional, Union
import pandas as pd
import logging
from pathlib import Path
import requests
import numpy as np
from functools import lru_cache

logger = logging.getLogger(__name__)

class GeneHarmonizer:
    """Class for harmonizing gene identifiers across datasets."""
    
    BIOMART_URL = "http://www.ensembl.org/biomart/martservice"
    
    def __init__(self, cache_dir: Optional[str] = None):
        """
        Initialize the gene harmonizer.
        
        Args:
            cache_dir: Directory to cache gene ID mappings
        """
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    @lru_cache(maxsize=1)
    def _get_id_mapping(self, source_type: str, target_type: str) -> Dict[str, str]:
        """
        Get gene ID mapping from biomart or cache.
        
        Args:
            source_type: Source ID type (e.g., 'ensembl_gene_id', 'hgnc_symbol')
            target_type: Target ID type
            
        Returns:
            Dict mapping source IDs to target IDs
        """
        cache_file = None
        if self.cache_dir:
            cache_file = self.cache_dir / f"{source_type}_to_{target_type}.parquet"
            if cache_file.exists():
                return pd.read_parquet(cache_file).set_index('source')['target'].to_dict()
        
        # Construct biomart query
        query = f"""<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="">
            <Dataset name="hsapiens_gene_ensembl" interface="default">
                <Attribute name="{source_type}"/>
                <Attribute name="{target_type}"/>
            </Dataset>
        </Query>
        """
        
        try:
            response = requests.get(self.BIOMART_URL, params={'query': query})
            response.raise_for_status()
            
            # Parse results
            mappings = {}
            for line in response.text.strip().split('\n'):
                source, target = line.split('\t')
                if source and target:  # Only keep non-empty mappings
                    mappings[source] = target
            
            # Cache results if cache_dir provided
            if cache_file:
                pd.DataFrame({'source': list(mappings.keys()),
                            'target': list(mappings.values())}).to_parquet(cache_file)
            
            return mappings
            
        except Exception as e:
            logger.error(f"Error fetching gene ID mappings: {str(e)}")
            raise
    
    def find_common_genes(self, datasets: List[pd.DataFrame],
                         id_types: List[str]) -> Set[str]:
        """
        Find genes common to all datasets.
        
        Args:
            datasets: List of gene expression matrices
            id_types: List of ID types for each dataset
            
        Returns:
            Set of common gene IDs in target type
        """
        common_genes = None
        
        for data, id_type in zip(datasets, id_types):
            # Get gene IDs from current dataset
            current_genes = set(data.columns if isinstance(data, pd.DataFrame)
                              else data.var_names)
            
            # Convert to target type if needed
            if id_type != 'ensembl_gene_id':
                mapping = self._get_id_mapping(id_type, 'ensembl_gene_id')
                current_genes = {mapping[g] for g in current_genes if g in mapping}
            
            # Update common genes
            if common_genes is None:
                common_genes = current_genes
            else:
                common_genes &= current_genes
        
        return common_genes if common_genes else set()
    
    def harmonize_matrices(self, datasets: List[Union[pd.DataFrame, pd.Series]],
                          id_types: List[str],
                          target_type: str = 'ensembl_gene_id') -> List[pd.DataFrame]:
        """
        Harmonize gene IDs across expression matrices and subset to common genes.
        
        Args:
            datasets: List of gene expression matrices
            id_types: List of ID types for each dataset
            target_type: Target ID type for harmonization
            
        Returns:
            List of harmonized matrices with common genes
        """
        # First find common genes
        common_genes = self.find_common_genes(datasets, id_types)
        logger.info(f"Found {len(common_genes)} common genes")
        
        harmonized_datasets = []
        for data, id_type in zip(datasets, id_types):
            # Convert to DataFrame if Series
            if isinstance(data, pd.Series):
                data = data.to_frame()
            
            # Convert current IDs to target type
            if id_type != target_type:
                mapping = self._get_id_mapping(id_type, target_type)
                
                # Create reverse mapping for common genes only
                rev_mapping = {v: k for k, v in mapping.items() 
                             if v in common_genes}
                
                # Subset and rename columns
                current_genes = [g for g in data.columns if g in mapping 
                               and mapping[g] in common_genes]
                harmonized = data[current_genes].copy()
                harmonized.columns = [mapping[g] for g in current_genes]
            else:
                # Just subset to common genes
                harmonized = data[list(common_genes)].copy()
            
            harmonized_datasets.append(harmonized)
        
        return harmonized_datasets
    
    def get_gene_info(self, gene_ids: List[str], 
                     id_type: str = 'ensembl_gene_id') -> pd.DataFrame:
        """
        Get additional information about genes.
        
        Args:
            gene_ids: List of gene IDs
            id_type: Type of input IDs
            
        Returns:
            DataFrame with gene information
        """
        # Construct biomart query for gene info
        id_list = ','.join(gene_ids)
        query = f"""<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="">
            <Dataset name="hsapiens_gene_ensembl" interface="default">
                <Filter name="{id_type}" value="{id_list}"/>
                <Attribute name="{id_type}"/>
                <Attribute name="external_gene_name"/>
                <Attribute name="chromosome_name"/>
                <Attribute name="start_position"/>
                <Attribute name="end_position"/>
                <Attribute name="strand"/>
                <Attribute name="gene_biotype"/>
                <Attribute name="description"/>
            </Dataset>
        </Query>
        """
        
        try:
            response = requests.get(self.BIOMART_URL, params={'query': query})
            response.raise_for_status()
            
            # Parse results into DataFrame
            info_df = pd.read_csv(pd.StringIO(response.text), sep='\t')
            return info_df
            
        except Exception as e:
            logger.error(f"Error fetching gene information: {str(e)}")
            raise 