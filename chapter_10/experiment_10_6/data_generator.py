#!/usr/bin/env python3
"""
Synthetic Epigenomic Data Generator
Generates realistic epigenomic track data for HMM segmentation testing
"""

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from scipy import stats

def generate_epigenomic_data(n_positions=10000, n_states=3, min_segment_length=50, 
                           chromosome="chr1", start_position=1000000, 
                           noise_level=0.1, seed=42):
    """
    Generate synthetic epigenomic data with defined chromatin states
    
    Parameters:
    - n_positions: Number of genomic positions
    - n_states: Number of chromatin states
    - min_segment_length: Minimum length of each segment
    - chromosome: Chromosome name
    - start_position: Starting genomic coordinate
    - noise_level: Amount of noise to add
    - seed: Random seed for reproducibility
    """
    
    np.random.seed(seed)
    
    # Define state characteristics (realistic epigenomic patterns)
    state_params = {
        0: {'name': 'Heterochromatin', 'mean': 0.5, 'std': 0.3, 'color': 'red'},
        1: {'name': 'Euchromatin', 'mean': 2.0, 'std': 0.5, 'color': 'green'},
        2: {'name': 'Active_Promoter', 'mean': 4.5, 'std': 0.8, 'color': 'blue'}
    }
    
    if n_states > 3:
        # Add more states for larger models
        for i in range(3, n_states):
            state_params[i] = {
                'name': f'State_{i}',
                'mean': 1.0 + i * 1.2,
                'std': 0.4 + i * 0.1,
                'color': plt.cm.tab10(i)
            }
    
    # Generate transition matrix (favor self-transitions)
    transition_matrix = np.full((n_states, n_states), 0.05)
    np.fill_diagonal(transition_matrix, 0.85)
    
    # Normalize rows
    transition_matrix = transition_matrix / transition_matrix.sum(axis=1, keepdims=True)
    
    # Generate state sequence using HMM
    states = []
    current_state = np.random.choice(n_states)
    segment_length = 0
    
    for i in range(n_positions):
        # Force minimum segment length
        if segment_length < min_segment_length:
            states.append(current_state)
            segment_length += 1
        else:
            # Allow state transition
            if np.random.random() < 0.1:  # 10% chance of transition
                current_state = np.random.choice(n_states, p=transition_matrix[current_state])
                segment_length = 0
            states.append(current_state)
            segment_length += 1
    
    # Generate observations based on states
    observations = []
    for state in states:
        params = state_params[state]
        # Generate from normal distribution with state-specific parameters
        obs = np.random.normal(params['mean'], params['std'])
        # Add some noise
        obs += np.random.normal(0, noise_level)
        # Ensure non-negative values (common in epigenomic data)
        obs = max(0, obs)
        observations.append(obs)
    
    # Create genomic positions
    positions = list(range(start_position, start_position + n_positions * 1000, 1000))
    
    # Create DataFrame
    df = pd.DataFrame({
        'chromosome': [chromosome] * n_positions,
        'position': positions,
        'signal': observations,
        'true_state': states
    })
    
    return df, state_params, transition_matrix

def add_realistic_features(df, noise_types=['outliers', 'missing', 'batch_effects']):
    """Add realistic features to make data more challenging"""
    
    df_modified = df.copy()
    n = len(df_modified)
    
    if 'outliers' in noise_types:
        # Add some outliers (5% of data)
        outlier_indices = np.random.choice(n, size=int(0.05 * n), replace=False)
        for idx in outlier_indices:
            df_modified.loc[idx, 'signal'] *= np.random.uniform(3, 8)
    
    if 'missing' in noise_types:
        # Add missing values (2% of data)
        missing_indices = np.random.choice(n, size=int(0.02 * n), replace=False)
        df_modified.loc[missing_indices, 'signal'] = np.nan
    
    if 'batch_effects' in noise_types:
        # Add batch effects (systematic shifts in different regions)
        batch_size = n // 4
        batches = [0, 1, 2, 3]
        batch_effects = [0, 0.2, -0.1, 0.3]
        
        for i, effect in enumerate(batch_effects):
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, n)
            df_modified.loc[start_idx:end_idx, 'signal'] += effect
    
    return df_modified

def create_visualization(df, state_params, output_file='epigenomic_data_visualization.png'):
    """Create visualization of the generated data"""
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 10))
    
    # Plot 1: Signal track
    ax1.plot(df['position'], df['signal'], alpha=0.7, linewidth=0.5)
    ax1.set_ylabel('Signal Intensity')
    ax1.set_title('Synthetic Epigenomic Signal')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: True states
    colors = [state_params[state]['color'] for state in df['true_state']]
    ax2.scatter(df['position'], df['true_state'], c=colors, alpha=0.6, s=1)
    ax2.set_ylabel('Chromatin State')
    ax2.set_title('True Chromatin States')
    ax2.set_yticks(range(len(state_params)))
    ax2.set_yticklabels([state_params[i]['name'] for i in range(len(state_params))])
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Signal colored by state
    for state in state_params:
        mask = df['true_state'] == state
        ax3.scatter(df.loc[mask, 'position'], df.loc[mask, 'signal'], 
                   c=state_params[state]['color'], alpha=0.6, s=1,
                   label=state_params[state]['name'])
    
    ax3.set_xlabel('Genomic Position')
    ax3.set_ylabel('Signal Intensity')
    ax3.set_title('Signal Intensity by Chromatin State')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Visualization saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate synthetic epigenomic data')
    parser.add_argument('--n-positions', type=int, default=5000,
                       help='Number of genomic positions (default: 5000)')
    parser.add_argument('--n-states', type=int, default=3,
                       help='Number of chromatin states (default: 3)')
    parser.add_argument('--min-segment-length', type=int, default=50,
                       help='Minimum segment length (default: 50)')
    parser.add_argument('--chromosome', default='chr1',
                       help='Chromosome name (default: chr1)')
    parser.add_argument('--start-position', type=int, default=1000000,
                       help='Starting genomic position (default: 1000000)')
    parser.add_argument('--noise-level', type=float, default=0.1,
                       help='Noise level (default: 0.1)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed (default: 42)')
    parser.add_argument('--output', default='epigenomic_data.csv',
                       help='Output CSV file (default: epigenomic_data.csv)')
    parser.add_argument('--add-noise', action='store_true',
                       help='Add realistic noise and artifacts')
    parser.add_argument('--visualize', action='store_true',
                       help='Create visualization plots')
    
    args = parser.parse_args()
    
    print("Generating synthetic epigenomic data...")
    print(f"Parameters: {args.n_positions} positions, {args.n_states} states")
    
    # Generate data
    df, state_params, transition_matrix = generate_epigenomic_data(
        n_positions=args.n_positions,
        n_states=args.n_states,
        min_segment_length=args.min_segment_length,
        chromosome=args.chromosome,
        start_position=args.start_position,
        noise_level=args.noise_level,
        seed=args.seed
    )
    
    # Add realistic noise if requested
    if args.add_noise:
        print("Adding realistic noise and artifacts...")
        df = add_realistic_features(df)
    
    # Remove rows with missing values for CSV output
    df_clean = df.dropna()
    
    # Save to CSV (without true_state column for the HMM tool)
    output_df = df_clean[['chromosome', 'position', 'signal']].copy()
    output_df.to_csv(args.output, index=False)
    
    # Save complete data with true states for evaluation
    evaluation_file = args.output.replace('.csv', '_with_true_states.csv')
    df_clean.to_csv(evaluation_file, index=False)
    
    print(f"Data saved to {args.output}")
    print(f"Data with true states saved to {evaluation_file}")
    
    # Create visualization if requested
    if args.visualize:
        create_visualization(df_clean, state_params)
    
    # Print summary statistics
    print("\n=== Data Summary ===")
    print(f"Total positions: {len(df_clean)}")
    print(f"Signal range: {df_clean['signal'].min():.3f} - {df_clean['signal'].max():.3f}")
    print(f"Signal mean: {df_clean['signal'].mean():.3f}")
    print(f"Signal std: {df_clean['signal'].std():.3f}")
    
    print("\nState distribution:")
    state_counts = df_clean['true_state'].value_counts().sort_index()
    for state, count in state_counts.items():
        percentage = count / len(df_clean) * 100
        print(f"  State {state} ({state_params[state]['name']}): {count} ({percentage:.1f}%)")
    
    print("\nTransition matrix:")
    print(transition_matrix)

if __name__ == "__main__":
    main()