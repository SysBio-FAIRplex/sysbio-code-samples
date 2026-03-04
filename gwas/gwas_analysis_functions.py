#!/usr/bin/env python3
"""
GWAS Analysis Functions for AD vs PD Comparison

This module contains functions for:
- Sample extraction and processing
- Quality assessment and lambda calculation
- Result filtering and analysis
- Statistical calculations

Author: [Your Name]
Date: [Date]
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import sys

def extract_samples():
    """
    Extract AD and PD samples from metadata files.
    
    Creates sample lists for cases and controls from both datasets.
    """
    print("Processing PD samples...")
    
    # Extract PD samples
    amp_pd_meta_df = pd.read_csv("/path/to/AMP_PD/metadata/AMPPD_v3_COV_wPHENOS_wGENOTOOLS.csv")
    amp_pd_control_subset_df = amp_pd_meta_df[amp_pd_meta_df['LATEST_DX'] == 'No PD Nor Other Neurological Disorder']
    amp_pd_case_subset_df = amp_pd_meta_df[amp_pd_meta_df['LATEST_DX'].isin(['Parkinson\'s Disease', 'Idiopathic PD'])]

    # Save PD sample lists
    amp_pd_control_subset_df[["ID","ID"]].to_csv("amp_pd_control_list", sep="\t", index=False)
    amp_pd_case_subset_df[["ID","ID"]].to_csv("amp_pd_case_list", sep="\t", index=False)

    # Extract AD samples
    print("Processing AD samples...")
    amp_ad_meta_df = pd.read_csv("/path/to/AMP_AD/metadata/threeCohorts.csv")
    amp_ad_meta_df['ID'] = amp_ad_meta_df['wgs_id'].astype(str).str.split('.').str[0]
    amp_ad_control_subset_df = amp_ad_meta_df[amp_ad_meta_df['clinAD'] == 0]
    amp_ad_case_subset_df = amp_ad_meta_df[amp_ad_meta_df['clinAD'] == 1]

    # Save AD sample lists
    amp_ad_control_subset_df[["ID","ID"]].to_csv("amp_ad_control_list", sep="\t", index=False)
    amp_ad_case_subset_df[["ID","ID"]].to_csv("amp_ad_case_list", sep="\t", index=False)

    print("Sample lists created successfully")
    
    return {
        'pd_controls': len(amp_pd_control_subset_df),
        'pd_cases': len(amp_pd_case_subset_df),
        'ad_controls': len(amp_ad_control_subset_df),
        'ad_cases': len(amp_ad_case_subset_df)
    }

def calculate_lambda(gwas_file, test_type="ADD", n_cases=None, n_controls=None):
    """
    Calculate genomic inflation factor (lambda) for GWAS results.
    
    Parameters:
    -----------
    gwas_file : str
        Path to GWAS results file (.assoc.logistic)
    test_type : str, default "ADD"
        Type of test to analyze (usually "ADD" for additive model)
    n_cases : int, optional
        Number of cases for scaled lambda calculation
    n_controls : int, optional
        Number of controls for scaled lambda calculation
    
    Returns:
    --------
    dict
        Dictionary containing lambda statistics
    """
    print(f"Analyzing GWAS results from: {gwas_file}")
    
    # Read GWAS results
    gwas_results = pd.read_csv(gwas_file, delim_whitespace=True, engine='c')
    gwas_results = gwas_results[gwas_results['TEST'] == test_type]
    
    # Convert P-values to chi-square statistics
    gwas_results['CHISQ'] = stats.chi2.ppf(1 - gwas_results['P'], 1)
    
    # Calculate lambda
    chi2 = gwas_results['CHISQ'].dropna()
    lambda_gc = np.median(chi2) / 0.455  # 0.455 is the median of chi2 distribution with 1 df
    
    # Calculate scaled lambda if sample sizes provided
    lambda_scaled = None
    if n_cases and n_controls:
        lambda_scaled = 1 + (lambda_gc - 1) * ((1/n_cases + 1/n_controls) / (1/1000 + 1/1000))
    
    # Extract significant results (P < 5e-8)
    significant_results = gwas_results[gwas_results['P'] < 0.00000005]
    
    return {
        'lambda_gc': lambda_gc,
        'lambda_scaled': lambda_scaled,
        'n_variants': len(gwas_results),
        'n_significant': len(significant_results),
        'significant_results': significant_results
    }

def assess_quality(cases_file, controls_file):
    """
    Perform comprehensive quality assessment of GWAS results.
    
    Parameters:
    -----------
    cases_file : str
        Path to cases comparison GWAS results
    controls_file : str
        Path to controls comparison GWAS results
    
    Returns:
    --------
    dict
        Quality assessment results
    """
    print("Performing quality assessment...")
    
    # Sample sizes (from the original analysis)
    cases_n_cases = 561
    cases_n_controls = 2542
    controls_n_cases = 832
    controls_n_controls = 3031
    
    # Analyze cases comparison
    print("\n=== Cases Comparison (AD cases vs PD cases) ===")
    cases_results = calculate_lambda(cases_file, n_cases=cases_n_cases, n_controls=cases_n_controls)
    
    print(f"Lambda: {cases_results['lambda_gc']:.4f}")
    print(f"Scaled Lambda: {cases_results['lambda_scaled']:.4f}")
    print(f"Number of significant variants: {cases_results['n_significant']}")
    
    # Save significant results
    if len(cases_results['significant_results']) > 0:
        cases_results['significant_results'].to_csv('merged_cases_AD_pheno_GWAS_significant_results.txt', 
                                                   sep='\t', index=False)
    
    # Analyze controls comparison
    print("\n=== Controls Comparison (AD controls vs PD controls) ===")
    controls_results = calculate_lambda(controls_file, n_cases=controls_n_cases, n_controls=controls_n_controls)
    
    print(f"Lambda: {controls_results['lambda_gc']:.4f}")
    print(f"Scaled Lambda: {controls_results['lambda_scaled']:.4f}")
    print(f"Number of significant variants: {controls_results['n_significant']}")
    
    # Save significant results
    if len(controls_results['significant_results']) > 0:
        controls_results['significant_results'].to_csv('merged_controls_AD_pheno_GWAS_significant_results.txt', 
                                                      sep='\t', index=False)
    
    print("\nQuality assessment complete!")
    
    return {
        'cases': cases_results,
        'controls': controls_results
    }

def filter_results(cases_file, controls_file, output_file='case_versus_case_filtered_results.csv'):
    """
    Filter case results to remove variants significantly associated in controls.
    
    Parameters:
    -----------
    cases_file : str
        Path to cases comparison GWAS results
    controls_file : str
        Path to controls comparison GWAS results
    output_file : str
        Output file for filtered results
    
    Returns:
    --------
    int
        Number of variants remaining after filtering
    """
    print("Filtering case results to remove control-associated variants...")
    
    # Read GWAS results
    cases_df = pd.read_csv(cases_file, delim_whitespace=True, engine='c')
    controls_df = pd.read_csv(controls_file, delim_whitespace=True, engine='c')
    
    # Filter for additive model results only
    filtered_cases = cases_df[cases_df['TEST'].isin(["ADD"])]
    filtered_controls = controls_df[controls_df['TEST'].isin(["ADD"])]
    
    # Get SNPs from controls where P > 1E-5 (not significantly associated in controls)
    insignificant_control_snps = filtered_controls[filtered_controls['P'] > 1E-5]['SNP']
    
    # Filter cases for those SNPs (remove variants significantly associated in controls)
    final_df = filtered_cases[filtered_cases['SNP'].isin(insignificant_control_snps)]
    
    # Export filtered results
    final_df.to_csv(output_file, index=False)
    print(f"Filtered results saved to {output_file}")
    print(f"Remaining variants: {len(final_df)}")
    
    return len(final_df)

def create_qq_plot(gwas_file, output_file=None, title="Q-Q Plot"):
    """
    Create a Q-Q plot for GWAS results.
    
    Parameters:
    -----------
    gwas_file : str
        Path to GWAS results file
    output_file : str, optional
        Output file for the plot
    title : str
        Title for the plot
    
    Returns:
    --------
    matplotlib.figure.Figure
        The Q-Q plot figure
    """
    try:
        import matplotlib.pyplot as plt
        
        # Read GWAS results
        gwas_results = pd.read_csv(gwas_file, delim_whitespace=True, engine='c')
        gwas_results = gwas_results[gwas_results['TEST'] == "ADD"]
        
        # Remove NaN p-values
        p_values = gwas_results['P'].dropna()
        
        # Create Q-Q plot
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Expected p-values
        expected = -np.log10(np.linspace(0, 1, len(p_values) + 1)[1:])
        expected.sort()
        
        # Observed p-values
        observed = -np.log10(p_values.sort_values())
        
        # Plot
        ax.scatter(expected, observed, alpha=0.6, s=1)
        ax.plot([0, max(expected)], [0, max(expected)], 'r--', alpha=0.8)
        
        ax.set_xlabel('Expected -log10(P)')
        ax.set_ylabel('Observed -log10(P)')
        ax.set_title(title)
        
        # Add lambda annotation
        lambda_gc = calculate_lambda(gwas_file)['lambda_gc']
        ax.text(0.05, 0.95, f'λ = {lambda_gc:.3f}', 
                transform=ax.transAxes, fontsize=12,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Q-Q plot saved to {output_file}")
        
        return fig
        
    except ImportError:
        print("matplotlib not available. Skipping Q-Q plot creation.")
        return None

def summarize_results(cases_file, controls_file, filtered_file=None):
    """
    Create a summary of GWAS results.
    
    Parameters:
    -----------
    cases_file : str
        Path to cases comparison GWAS results
    controls_file : str
        Path to controls comparison GWAS results
    filtered_file : str, optional
        Path to filtered results file
    
    Returns:
    --------
    dict
        Summary statistics
    """
    print("Creating results summary...")
    
    # Read files
    cases_df = pd.read_csv(cases_file, delim_whitespace=True, engine='c')
    controls_df = pd.read_csv(controls_file, delim_whitespace=True, engine='c')
    
    # Filter for additive model
    cases_add = cases_df[cases_df['TEST'] == "ADD"]
    controls_add = controls_df[controls_df['TEST'] == "ADD"]
    
    # Calculate summary statistics
    summary = {
        'cases': {
            'total_variants': len(cases_add),
            'significant_variants': len(cases_add[cases_add['P'] < 5e-8]),
            'suggestive_variants': len(cases_add[cases_add['P'] < 1e-5]),
            'min_p_value': cases_add['P'].min(),
            'lambda': calculate_lambda(cases_file)['lambda_gc']
        },
        'controls': {
            'total_variants': len(controls_add),
            'significant_variants': len(controls_add[controls_add['P'] < 5e-8]),
            'suggestive_variants': len(controls_add[controls_add['P'] < 1e-5]),
            'min_p_value': controls_add['P'].min(),
            'lambda': calculate_lambda(controls_file)['lambda_gc']
        }
    }
    
    # Add filtered results if available
    if filtered_file and os.path.exists(filtered_file):
        filtered_df = pd.read_csv(filtered_file)
        summary['filtered'] = {
            'total_variants': len(filtered_df),
            'significant_variants': len(filtered_df[filtered_df['P'] < 5e-8]),
            'suggestive_variants': len(filtered_df[filtered_df['P'] < 1e-5]),
            'min_p_value': filtered_df['P'].min()
        }
    
    # Print summary
    print("\n=== RESULTS SUMMARY ===")
    print(f"\nCases Comparison:")
    print(f"  Total variants: {summary['cases']['total_variants']:,}")
    print(f"  Significant variants (P < 5e-8): {summary['cases']['significant_variants']}")
    print(f"  Suggestive variants (P < 1e-5): {summary['cases']['suggestive_variants']}")
    print(f"  Minimum P-value: {summary['cases']['min_p_value']:.2e}")
    print(f"  Lambda: {summary['cases']['lambda']:.4f}")
    
    print(f"\nControls Comparison:")
    print(f"  Total variants: {summary['controls']['total_variants']:,}")
    print(f"  Significant variants (P < 5e-8): {summary['controls']['significant_variants']}")
    print(f"  Suggestive variants (P < 1e-5): {summary['controls']['suggestive_variants']}")
    print(f"  Minimum P-value: {summary['controls']['min_p_value']:.2e}")
    print(f"  Lambda: {summary['controls']['lambda']:.4f}")
    
    if 'filtered' in summary:
        print(f"\nFiltered Results:")
        print(f"  Total variants: {summary['filtered']['total_variants']:,}")
        print(f"  Significant variants (P < 5e-8): {summary['filtered']['significant_variants']}")
        print(f"  Suggestive variants (P < 1e-5): {summary['filtered']['suggestive_variants']}")
        print(f"  Minimum P-value: {summary['filtered']['min_p_value']:.2e}")
    
    return summary

if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 1:
        command = sys.argv[1]
        
        if command == "extract_samples":
            extract_samples()
        elif command == "quality_assessment":
            if len(sys.argv) >= 4:
                assess_quality(sys.argv[2], sys.argv[3])
            else:
                print("Usage: python gwas_analysis_functions.py quality_assessment <cases_file> <controls_file>")
        elif command == "filter_results":
            if len(sys.argv) >= 4:
                filter_results(sys.argv[2], sys.argv[3])
            else:
                print("Usage: python gwas_analysis_functions.py filter_results <cases_file> <controls_file>")
        elif command == "summarize":
            if len(sys.argv) >= 4:
                summarize_results(sys.argv[2], sys.argv[3])
            else:
                print("Usage: python gwas_analysis_functions.py summarize <cases_file> <controls_file>")
        else:
            print("Unknown command. Available commands: extract_samples, quality_assessment, filter_results, summarize")
    else:
        print("GWAS Analysis Functions Module")
        print("Available functions:")
        print("- extract_samples(): Extract AD and PD samples")
        print("- assess_quality(cases_file, controls_file): Quality assessment")
        print("- filter_results(cases_file, controls_file): Filter results")
        print("- summarize_results(cases_file, controls_file): Create summary")
        print("- create_qq_plot(gwas_file): Create Q-Q plot") 