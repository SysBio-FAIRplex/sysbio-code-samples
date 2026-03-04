import pandas as pd
import numpy as np
import anndata as ad
from typing import Optional, Dict, List, Union
import logging
from pathlib import Path
import allel
import pysam

from .omics_processor import OmicsProcessor

logger = logging.getLogger(__name__)

class GenomicsProcessor(OmicsProcessor):
    """Processor for genomics data."""
    
    def __init__(self, data_dir: str, cache_dir: Optional[str] = None,
                 reference_genome: Optional[str] = None):
        """
        Initialize the genomics processor.
        
        Args:
            data_dir: Directory for data storage
            cache_dir: Directory for caching
            reference_genome: Path to reference genome FASTA file
        """
        super().__init__(data_dir, cache_dir)
        self.reference_genome = reference_genome
        self.fasta = pysam.FastaFile(reference_genome) if reference_genome else None
        
    def process_data(self, data: Union[pd.DataFrame, ad.AnnData],
                    processing_steps: List[Dict]) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Process genomics data according to specified steps.
        
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
                if step_type == 'filter_variants':
                    # Filter variants by quality and frequency
                    min_qual = params.get('min_quality', 20)
                    min_af = params.get('min_allele_freq', 0.01)
                    
                    if isinstance(processed, ad.AnnData):
                        qual_mask = processed.var['QUAL'] >= min_qual
                        af_mask = processed.var['AF'] >= min_af
                        processed = processed[:, qual_mask & af_mask]
                    else:
                        qual_mask = processed['QUAL'] >= min_qual
                        af_mask = processed['AF'] >= min_af
                        processed = processed[qual_mask & af_mask]
                        
                elif step_type == 'annotate_variants':
                    if self.reference_genome:
                        self._annotate_variants(processed)
                        
                elif step_type == 'calculate_frequencies':
                    if isinstance(processed, ad.AnnData):
                        gt_matrix = processed.X
                        af = np.mean(gt_matrix, axis=0)
                        processed.var['AF'] = af
                    else:
                        gt_cols = [col for col in processed.columns if col.startswith('GT')]
                        if gt_cols:
                            processed['AF'] = processed[gt_cols].mean(axis=1)
                            
            except Exception as e:
                logger.error(f"Error in processing step {step_name}: {str(e)}")
                raise
                
        return processed
        
    def calculate_quality_metrics(self, data: Union[pd.DataFrame, ad.AnnData]) -> Dict:
        """
        Calculate genomics-specific quality metrics.
        
        Args:
            data: Input data
            
        Returns:
            Dict with quality metrics
        """
        metrics = self.validate_data_quality(data, 'genomics')
        
        # Add genomics-specific metrics
        if isinstance(data, ad.AnnData):
            metrics.update({
                'total_variants': data.n_vars,
                'variants_per_chromosome': data.var['chromosome'].value_counts().to_dict(),
                'transition_transversion_ratio': self._calculate_ti_tv_ratio(data.var),
                'heterozygosity': np.mean(data.X == 1) if data.X.dtype in [np.int32, np.int64] else None
            })
        else:
            metrics.update({
                'total_variants': len(data),
                'variants_per_chromosome': data['chromosome'].value_counts().to_dict(),
                'transition_transversion_ratio': self._calculate_ti_tv_ratio(data),
                'heterozygosity': self._calculate_heterozygosity(data)
            })
            
        return metrics
        
    def _calculate_ti_tv_ratio(self, variants: pd.DataFrame) -> float:
        """Calculate transition/transversion ratio."""
        transitions = ['AG', 'GA', 'CT', 'TC']
        transversions = ['AC', 'CA', 'AT', 'TA', 'GC', 'CG', 'GT', 'TG']
        
        if 'reference' in variants.columns and 'alternate' in variants.columns:
            changes = variants['reference'] + variants['alternate']
            ti_count = sum(changes.str.contains('|'.join(transitions)))
            tv_count = sum(changes.str.contains('|'.join(transversions)))
            
            return ti_count / tv_count if tv_count > 0 else 0
        return 0
        
    def _calculate_heterozygosity(self, data: pd.DataFrame) -> float:
        """Calculate heterozygosity rate."""
        gt_cols = [col for col in data.columns if col.startswith('GT')]
        if gt_cols:
            genotypes = data[gt_cols]
            het_count = np.sum(genotypes == 1)
            total_count = np.prod(genotypes.shape)
            return het_count / total_count
        return 0
        
    def _annotate_variants(self, data: Union[pd.DataFrame, ad.AnnData]):
        """Add variant annotations using reference genome."""
        if not self.fasta:
            raise ValueError("Reference genome not provided")
            
        if isinstance(data, ad.AnnData):
            variants = data.var
        else:
            variants = data
            
        # Add sequence context
        context_size = 10  # bases on each side
        contexts = []
        
        for _, row in variants.iterrows():
            chrom = row['chromosome']
            pos = row['position']
            try:
                context = self.fasta.fetch(chrom,
                                         max(0, pos - context_size - 1),
                                         pos + context_size)
                contexts.append(context)
            except Exception:
                contexts.append('')
                
        if isinstance(data, ad.AnnData):
            data.var['sequence_context'] = contexts
        else:
            data['sequence_context'] = contexts
            
    def convert_coordinates(self, data: Union[pd.DataFrame, ad.AnnData],
                          source_build: str,
                          target_build: str) -> Union[pd.DataFrame, ad.AnnData]:
        """
        Convert genomic coordinates between genome builds.
        
        Args:
            data: Input data
            source_build: Source genome build
            target_build: Target genome build
            
        Returns:
            Data with converted coordinates
        """
        # TODO: Implement coordinate conversion using liftover
        raise NotImplementedError("Coordinate conversion not yet implemented")
        
    def merge_variants(self, datasets: List[Union[pd.DataFrame, ad.AnnData]],
                      merge_strategy: str = 'intersection') -> Union[pd.DataFrame, ad.AnnData]:
        """
        Merge variant datasets.
        
        Args:
            datasets: List of variant datasets
            merge_strategy: Strategy for merging ('intersection' or 'union')
            
        Returns:
            Merged dataset
        """
        # Convert all datasets to DataFrame format
        dfs = []
        for data in datasets:
            if isinstance(data, ad.AnnData):
                df = pd.DataFrame(data.X, columns=data.var_names, index=data.obs_names)
                df = pd.concat([df, data.var], axis=1)
                dfs.append(df)
            else:
                dfs.append(data)
                
        # Create variant keys
        for df in dfs:
            df['variant_key'] = df['chromosome'] + '_' + df['position'].astype(str) + '_' + \
                              df['reference'] + '_' + df['alternate']
                              
        # Merge based on strategy
        if merge_strategy == 'intersection':
            common_variants = set.intersection(*[set(df['variant_key']) for df in dfs])
            merged = pd.concat([df[df['variant_key'].isin(common_variants)] for df in dfs])
        else:  # union
            merged = pd.concat(dfs)
            
        return merged 