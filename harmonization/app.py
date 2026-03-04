import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path
import tempfile
import os

from src.data_processors.transcriptomics_processor import TranscriptomicsProcessor
from src.data_processors.proteomics_processor import ProteomicsProcessor
from src.data_processors.genomics_processor import GenomicsProcessor

# Set page config
st.set_page_config(
    page_title="Multi-Omics Data Processor",
    page_icon="🧬",
    layout="wide"
)

# Initialize session state
if 'processed_data' not in st.session_state:
    st.session_state.processed_data = None
if 'quality_metrics' not in st.session_state:
    st.session_state.quality_metrics = None

def setup_directories():
    """Create temporary directories for data processing."""
    temp_dir = tempfile.mkdtemp()
    data_dir = Path(temp_dir) / "data"
    cache_dir = Path(temp_dir) / "cache"
    
    data_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    return str(data_dir), str(cache_dir)

def process_transcriptomics(uploaded_files, processing_params):
    """Process transcriptomics data."""
    data_dir, cache_dir = setup_directories()
    processor = TranscriptomicsProcessor(data_dir=data_dir, cache_dir=cache_dir)
    
    # Save uploaded files
    datasets = []
    for file in uploaded_files:
        file_path = Path(data_dir) / file.name
        with open(file_path, 'wb') as f:
            f.write(file.getvalue())
        datasets.append(pd.read_csv(file_path))
    
    # Process each dataset
    processed_datasets = []
    metrics_list = []
    
    for data in datasets:
        processed = processor.process_data(data, [
            {
                'name': 'filter_genes',
                'type': 'filter_genes',
                'params': {'min_samples': processing_params['min_samples']}
            },
            {
                'name': 'normalize',
                'type': 'normalize',
                'params': {
                    'method': processing_params['norm_method'],
                    'scale_factor': processing_params['scale_factor']
                }
            }
        ])
        processed_datasets.append(processed)
        metrics_list.append(processor.calculate_quality_metrics(processed))
    
    # Harmonize if multiple datasets
    if len(processed_datasets) > 1:
        processed_datasets = processor.harmonize_datasets(
            processed_datasets,
            id_types=['ensembl_gene_id', 'symbol'],
            target_type='ensembl_gene_id'
        )
    
    return processed_datasets, metrics_list

def process_proteomics(uploaded_file, processing_params):
    """Process proteomics data."""
    data_dir, cache_dir = setup_directories()
    processor = ProteomicsProcessor(data_dir=data_dir, cache_dir=cache_dir)
    
    # Save and load data
    file_path = Path(data_dir) / uploaded_file.name
    with open(file_path, 'wb') as f:
        f.write(uploaded_file.getvalue())
    data = pd.read_csv(file_path)
    
    # Process data
    processed = processor.process_data(data, [
        {
            'name': 'filter_proteins',
            'type': 'filter_proteins',
            'params': {'min_samples': processing_params['min_samples']}
        },
        {
            'name': 'normalize',
            'type': 'normalize',
            'params': {'method': processing_params['norm_method']}
        }
    ])
    
    # Process modifications if applicable
    if 'modification_type' in data.columns:
        processed = processor.process_modifications(processed, 'phospho')
        mod_stats = processor.quantify_modifications(processed)
    else:
        mod_stats = None
    
    metrics = processor.calculate_quality_metrics(processed)
    
    return processed, metrics, mod_stats

def process_genomics(uploaded_file, processing_params):
    """Process genomics data."""
    data_dir, cache_dir = setup_directories()
    processor = GenomicsProcessor(
        data_dir=data_dir,
        cache_dir=cache_dir,
        reference_genome=None  # For demo, we'll skip reference genome
    )
    
    # Save and load data
    file_path = Path(data_dir) / uploaded_file.name
    with open(file_path, 'wb') as f:
        f.write(uploaded_file.getvalue())
    data = pd.read_csv(file_path)
    
    # Process data
    processed = processor.process_data(data, [
        {
            'name': 'filter_variants',
            'type': 'filter_variants',
            'params': {
                'min_quality': processing_params['min_quality'],
                'min_allele_freq': processing_params['min_af']
            }
        }
    ])
    
    metrics = processor.calculate_quality_metrics(processed)
    
    return processed, metrics

def plot_quality_metrics(metrics, data_type):
    """Create quality metric visualizations."""
    if data_type == 'transcriptomics':
        fig = px.bar(
            x=['Sample Count', 'Feature Count', 'Missing Values'],
            y=[metrics['sample_count'], metrics['feature_count'], metrics['missing_values']],
            title='Basic Quality Metrics'
        )
        st.plotly_chart(fig)
        
        # Plot mitochondrial content
        if 'percent_mito' in metrics:
            st.metric('Mitochondrial Content', f"{metrics['percent_mito']*100:.2f}%")
    
    elif data_type == 'proteomics':
        fig = px.bar(
            x=['Sample Count', 'Protein Count', 'Missing Rate (%)'],
            y=[metrics['sample_count'], metrics['feature_count'], 
               metrics['missing_value_rate']*100],
            title='Proteomics Quality Metrics'
        )
        st.plotly_chart(fig)
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric('Median Intensity', f"{metrics['median_intensity']:,.0f}")
        with col2:
            st.metric('Dynamic Range', f"{metrics['dynamic_range']:.2f}")
    
    elif data_type == 'genomics':
        # Plot variant distribution by chromosome
        chrom_data = pd.DataFrame(
            metrics['variants_per_chromosome'].items(),
            columns=['Chromosome', 'Variant Count']
        )
        fig = px.bar(chrom_data, x='Chromosome', y='Variant Count',
                    title='Variants per Chromosome')
        st.plotly_chart(fig)
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric('Ti/Tv Ratio', f"{metrics['transition_transversion_ratio']:.2f}")
        with col2:
            st.metric('Heterozygosity', f"{metrics['heterozygosity']:.3f}")

def main():
    st.title("🧬 Multi-Omics Data Processor")
    
    # Sidebar for data type selection and parameters
    st.sidebar.header("Settings")
    data_type = st.sidebar.selectbox(
        "Select Data Type",
        ["Transcriptomics", "Proteomics", "Genomics"]
    )
    
    # Data type specific parameters
    if data_type == "Transcriptomics":
        st.sidebar.subheader("Processing Parameters")
        processing_params = {
            'min_samples': st.sidebar.slider("Minimum Samples", 1, 10, 5),
            'norm_method': st.sidebar.selectbox(
                "Normalization Method",
                ["log1p", "log2", "zscore"]
            ),
            'scale_factor': st.sidebar.number_input(
                "Scale Factor",
                value=1e6,
                format="%e"
            )
        }
        
        uploaded_files = st.file_uploader(
            "Upload RNA-seq Data Files",
            accept_multiple_files=True,
            type=['csv', 'tsv']
        )
        
        if uploaded_files:
            if st.button("Process Data"):
                with st.spinner("Processing transcriptomics data..."):
                    processed_data, metrics = process_transcriptomics(
                        uploaded_files,
                        processing_params
                    )
                    st.session_state.processed_data = processed_data
                    st.session_state.quality_metrics = metrics
                    
                st.success("Processing complete!")
                
                # Display results
                for i, metrics in enumerate(st.session_state.quality_metrics):
                    st.subheader(f"Dataset {i+1} Results")
                    plot_quality_metrics(metrics, 'transcriptomics')
                    
                # Download processed data
                for i, data in enumerate(st.session_state.processed_data):
                    st.download_button(
                        f"Download Processed Dataset {i+1}",
                        data.to_csv(index=False).encode('utf-8'),
                        f"processed_transcriptomics_{i+1}.csv",
                        "text/csv"
                    )
    
    elif data_type == "Proteomics":
        st.sidebar.subheader("Processing Parameters")
        processing_params = {
            'min_samples': st.sidebar.slider("Minimum Samples", 1, 10, 3),
            'norm_method': st.sidebar.selectbox(
                "Normalization Method",
                ["log2", "zscore"]
            )
        }
        
        uploaded_file = st.file_uploader(
            "Upload Proteomics Data File",
            type=['csv', 'tsv']
        )
        
        if uploaded_file:
            if st.button("Process Data"):
                with st.spinner("Processing proteomics data..."):
                    processed_data, metrics, mod_stats = process_proteomics(
                        uploaded_file,
                        processing_params
                    )
                    st.session_state.processed_data = processed_data
                    st.session_state.quality_metrics = metrics
                    
                st.success("Processing complete!")
                
                # Display results
                st.subheader("Quality Metrics")
                plot_quality_metrics(metrics, 'proteomics')
                
                if mod_stats:
                    st.subheader("Modification Statistics")
                    st.json(mod_stats)
                
                # Download processed data
                st.download_button(
                    "Download Processed Data",
                    st.session_state.processed_data.to_csv(index=False).encode('utf-8'),
                    "processed_proteomics.csv",
                    "text/csv"
                )
    
    else:  # Genomics
        st.sidebar.subheader("Processing Parameters")
        processing_params = {
            'min_quality': st.sidebar.slider("Minimum Quality Score", 0, 100, 30),
            'min_af': st.sidebar.slider(
                "Minimum Allele Frequency",
                0.0, 1.0, 0.01,
                step=0.01
            )
        }
        
        uploaded_file = st.file_uploader(
            "Upload Genomics Data File",
            type=['csv', 'tsv']
        )
        
        if uploaded_file:
            if st.button("Process Data"):
                with st.spinner("Processing genomics data..."):
                    processed_data, metrics = process_genomics(
                        uploaded_file,
                        processing_params
                    )
                    st.session_state.processed_data = processed_data
                    st.session_state.quality_metrics = metrics
                    
                st.success("Processing complete!")
                
                # Display results
                st.subheader("Quality Metrics")
                plot_quality_metrics(metrics, 'genomics')
                
                # Download processed data
                st.download_button(
                    "Download Processed Data",
                    st.session_state.processed_data.to_csv(index=False).encode('utf-8'),
                    "processed_variants.csv",
                    "text/csv"
                )

if __name__ == "__main__":
    main() 