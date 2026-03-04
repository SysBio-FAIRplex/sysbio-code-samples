import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import logging
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import LabelEncoder
import re
from collections import Counter
import joblib
from pathlib import Path

logger = logging.getLogger(__name__)

class SchemaInferrer:
    """
    Intelligent schema inference for multi-omics biomedical data using ML techniques.
    """
    
    # Common patterns in biomedical data
    PATTERNS = {
        # Gene/Transcript patterns
        'ensembl_gene': r'^ENS[GTP]\d{11}$',
        'gene_symbol': r'^[A-Za-z0-9-]+$',
        'refseq_gene': r'^[NX][MR]_\d+$',
        
        # Protein patterns
        'uniprot_id': r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$',
        'protein_symbol': r'^[A-Z][A-Za-z0-9]*$',
        'peptide_sequence': r'^[ACDEFGHIKLMNPQRSTVWY]+$',
        
        # Genomic patterns
        'chromosome': r'^chr[0-9XYM]{1,2}$',
        'genomic_position': r'^\d+$',
        'variant_id': r'^rs\d+$',
        'allele': r'^[ACGT]+$',
        
        # Common data types
        'p_value': r'^[0-9.e-]+$',
        'sample_id': r'^[A-Za-z0-9-_]+$',
        'date': r'^\d{4}-\d{2}-\d{2}$'
    }
    
    # Common column name variations by data type
    COLUMN_VARIATIONS = {
        # Subject/Sample information
        'subject_id': ['subject', 'participant', 'individual', 'patient', 'donor'],
        'sample_id': ['sample', 'specimen', 'biospecimen', 'aliquot'],
        'condition': ['disease', 'phenotype', 'status', 'diagnosis'],
        'age': ['age_at', 'age_years', 'patient_age'],
        'sex': ['gender', 'biological_sex'],
        'tissue': ['tissue_type', 'sample_site', 'biopsy_site'],
        
        # Transcriptomics
        'gene_id': ['gene', 'ensembl', 'symbol', 'gene_name'],
        'expression': ['expr', 'counts', 'tpm', 'fpkm', 'abundance'],
        'transcript': ['transcript_id', 'isoform', 'splice_variant'],
        
        # Proteomics
        'protein_id': ['protein', 'uniprot', 'protein_name', 'accession'],
        'peptide': ['peptide_seq', 'amino_acid_sequence', 'aa_sequence'],
        'intensity': ['abundance', 'concentration', 'quantity', 'intensity'],
        'modification': ['ptm', 'phospho', 'ubiq', 'acetyl'],
        
        # Genomics
        'variant': ['snp', 'indel', 'mutation', 'variant_id'],
        'position': ['pos', 'coordinate', 'location', 'genomic_pos'],
        'reference': ['ref', 'reference_allele'],
        'alternate': ['alt', 'alternate_allele'],
        'quality': ['qual', 'quality_score', 'confidence']
    }
    
    # Data type specific metadata
    METADATA_TYPES = {
        'transcriptomics': {
            'required_columns': ['gene_id'],
            'value_type': 'expression',
            'common_transformations': ['log2', 'log10', 'zscore']
        },
        'proteomics': {
            'required_columns': ['protein_id'],
            'value_type': 'intensity',
            'common_transformations': ['log2', 'log10', 'zscore']
        },
        'genomics': {
            'required_columns': ['chromosome', 'position', 'reference', 'alternate'],
            'value_type': 'categorical',
            'common_transformations': []
        }
    }
    
    def __init__(self, model_dir: Optional[str] = None):
        """
        Initialize the schema inferrer.
        
        Args:
            model_dir: Directory to save/load trained models
        """
        self.model_dir = Path(model_dir) if model_dir else None
        if self.model_dir:
            self.model_dir.mkdir(parents=True, exist_ok=True)
            
        self.column_classifier = None
        self.tfidf = None
        self.label_encoder = None
        self.trained = False
        
    def _extract_column_features(self, column_name: str) -> List[str]:
        """Extract features from column name."""
        features = []
        
        # Split on common separators
        parts = re.split(r'[_\s-]', column_name.lower())
        features.extend(parts)
        
        # Check for common patterns
        for pattern_name, pattern in self.PATTERNS.items():
            if re.match(pattern, column_name):
                features.append(f'pattern_{pattern_name}')
                
        # Check for common variations
        for standard_name, variations in self.COLUMN_VARIATIONS.items():
            if any(var in column_name.lower() for var in variations):
                features.append(f'variation_{standard_name}')
                
        return features
        
    def _infer_data_type(self, series: pd.Series) -> str:
        """Infer the data type of a column."""
        # Get sample of non-null values
        sample = series.dropna().head(100)
        if len(sample) == 0:
            return 'unknown'
            
        # Check if all values match any of our patterns
        for pattern_name, pattern in self.PATTERNS.items():
            if all(isinstance(x, str) and re.match(pattern, str(x)) for x in sample):
                return pattern_name
                
        # Check numeric types
        if pd.api.types.is_numeric_dtype(series):
            if all(float(x).is_integer() for x in sample):
                return 'integer'
            return 'float'
            
        # Check categorical
        if len(series.unique()) / len(series) < 0.1:
            return 'categorical'
            
        return 'string'
        
    def train(self, training_data: List[Tuple[str, str]]):
        """
        Train the column classifier on labeled examples.
        
        Args:
            training_data: List of (column_name, column_type) tuples
        """
        # Extract features from column names
        column_features = [self._extract_column_features(name) for name, _ in training_data]
        
        # Convert to TF-IDF features
        self.tfidf = TfidfVectorizer()
        X = self.tfidf.fit_transform([' '.join(features) for features in column_features])
        
        # Encode labels
        self.label_encoder = LabelEncoder()
        y = self.label_encoder.fit_transform([col_type for _, col_type in training_data])
        
        # Train classifier
        self.column_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
        self.column_classifier.fit(X, y)
        self.trained = True
        
        # Save models if directory provided
        if self.model_dir:
            joblib.dump(self.column_classifier, self.model_dir / 'column_classifier.joblib')
            joblib.dump(self.tfidf, self.model_dir / 'tfidf.joblib')
            joblib.dump(self.label_encoder, self.model_dir / 'label_encoder.joblib')
            
    def load_models(self):
        """Load trained models from disk."""
        if not self.model_dir:
            raise ValueError("No model directory specified")
            
        self.column_classifier = joblib.load(self.model_dir / 'column_classifier.joblib')
        self.tfidf = joblib.load(self.model_dir / 'tfidf.joblib')
        self.label_encoder = joblib.load(self.model_dir / 'label_encoder.joblib')
        self.trained = True
        
    def infer_schema(self, df: pd.DataFrame) -> Dict:
        """
        Infer the schema of a DataFrame.
        
        Args:
            df: Input DataFrame
            
        Returns:
            Dict containing inferred column types and metadata
        """
        schema = {
            'columns': {},
            'metadata': {
                'num_rows': len(df),
                'num_cols': len(df.columns),
                'missing_values': df.isnull().sum().to_dict()
            }
        }
        
        for col in df.columns:
            col_features = self._extract_column_features(col)
            
            # Use trained classifier if available
            if self.trained:
                X = self.tfidf.transform([' '.join(col_features)])
                predicted_type = self.label_encoder.inverse_transform(
                    self.column_classifier.predict(X)
                )[0]
            else:
                # Use heuristic approach
                predicted_type = 'unknown'
                for standard_name, variations in self.COLUMN_VARIATIONS.items():
                    if any(var in col.lower() for var in variations):
                        predicted_type = standard_name
                        break
                        
            # Get data type
            data_type = self._infer_data_type(df[col])
            
            schema['columns'][col] = {
                'predicted_type': predicted_type,
                'data_type': data_type,
                'unique_values': len(df[col].unique()),
                'example_values': df[col].dropna().head(5).tolist()
            }
            
        return schema
        
    def infer_omics_type(self, df: pd.DataFrame) -> str:
        """
        Infer the type of omics data based on column patterns and content.
        
        Args:
            df: Input DataFrame
            
        Returns:
            String indicating the omics type ('transcriptomics', 'proteomics', 'genomics', or 'unknown')
        """
        # Get column types
        schema = self.infer_schema(df)
        column_types = {col: info['predicted_type'] for col, info in schema['columns'].items()}
        
        # Check for characteristic patterns of each omics type
        if any(t == 'gene_id' for t in column_types.values()) and \
           any(t == 'expression' for t in column_types.values()):
            return 'transcriptomics'
        
        if any(t in ['protein_id', 'peptide'] for t in column_types.values()) and \
           any(t == 'intensity' for t in column_types.values()):
            return 'proteomics'
        
        if all(req in column_types.values() for req in ['chromosome', 'position', 'reference', 'alternate']):
            return 'genomics'
        
        return 'unknown'
        
    def suggest_harmonization(self, schemas: List[Dict], omics_type: Optional[str] = None) -> Dict:
        """
        Suggest harmonization steps for multiple datasets.
        
        Args:
            schemas: List of schema dictionaries from infer_schema
            omics_type: Optional type of omics data to guide harmonization
            
        Returns:
            Dict with harmonization suggestions
        """
        suggestions = {
            'common_columns': [],
            'similar_columns': [],
            'type_conflicts': [],
            'harmonization_steps': [],
            'data_type_specific': {}
        }
        
        # If omics_type is provided, add specific requirements
        if omics_type and omics_type in self.METADATA_TYPES:
            metadata = self.METADATA_TYPES[omics_type]
            suggestions['data_type_specific'] = {
                'required_columns': metadata['required_columns'],
                'value_type': metadata['value_type'],
                'suggested_transformations': metadata['common_transformations']
            }
        
        # Find common columns across all schemas
        all_columns = [set(schema['columns'].keys()) for schema in schemas]
        common_columns = set.intersection(*all_columns)
        suggestions['common_columns'] = list(common_columns)
        
        # Find similar columns that might need harmonization
        column_groups = {}
        for schema in schemas:
            for col, info in schema['columns'].items():
                if col not in common_columns:
                    pred_type = info['predicted_type']
                    if pred_type not in column_groups:
                        column_groups[pred_type] = set()
                    column_groups[pred_type].add(col)
                    
        suggestions['similar_columns'] = [
            list(cols) for cols in column_groups.values() if len(cols) > 1
        ]
        
        # Check for type conflicts in common columns
        for col in common_columns:
            types = [schema['columns'][col]['data_type'] for schema in schemas]
            if len(set(types)) > 1:
                suggestions['type_conflicts'].append({
                    'column': col,
                    'types': types
                })
                
        # Generate harmonization steps
        for conflict in suggestions['type_conflicts']:
            col = conflict['column']
            types = conflict['types']
            
            # Suggest appropriate conversion
            if 'float' in types and 'integer' in types:
                suggestions['harmonization_steps'].append({
                    'column': col,
                    'action': 'convert_to_float',
                    'reason': 'Mixed numeric types'
                })
            elif 'categorical' in types:
                suggestions['harmonization_steps'].append({
                    'column': col,
                    'action': 'standardize_categories',
                    'reason': 'Inconsistent categorical values'
                })
                
        return suggestions
        
    def validate_requirements(self, df: pd.DataFrame, omics_type: str) -> Dict:
        """
        Validate that the dataset meets requirements for the specified omics type.
        
        Args:
            df: Input DataFrame
            omics_type: Type of omics data
            
        Returns:
            Dict with validation results
        """
        if omics_type not in self.METADATA_TYPES:
            raise ValueError(f"Unknown omics type: {omics_type}")
            
        metadata = self.METADATA_TYPES[omics_type]
        schema = self.infer_schema(df)
        
        validation = {
            'is_valid': True,
            'missing_required': [],
            'warnings': []
        }
        
        # Check required columns
        column_types = {col: info['predicted_type'] for col, info in schema['columns'].items()}
        for required in metadata['required_columns']:
            if required not in column_types.values():
                validation['is_valid'] = False
                validation['missing_required'].append(required)
        
        # Check value type
        value_cols = [col for col, info in schema['columns'].items() 
                     if info['predicted_type'] == metadata['value_type']]
        if not value_cols:
            validation['warnings'].append(f"No columns of type {metadata['value_type']} found")
        
        return validation 