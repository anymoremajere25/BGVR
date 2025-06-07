#!/usr/bin/env python3
"""
Evaluation script for HMM segmentation results
Compares predicted states with ground truth and generates metrics
"""

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, confusion_matrix
from sklearn.metrics import accuracy_score, precision_recall_fscore_support
import argparse
from scipy.optimize import linear_sum_assignment
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

def load_segmentation_results(results_file):
    """Load segmentation results from JSON file"""
    with open(results_file, 'r') as f:
        results = json.load(f)
    return results

def load_truth_data(truth_file):
    """Load ground truth data"""
    df = pd.read_csv(truth_file)
    return df

def align_states(true_states, pred_states):
    """
    Align predicted states with true states using Hungarian algorithm
    to account for label switching problem in unsupervised learning
    """
    # Create confusion matrix
    n_true = len(set(true_states))
    n_pred = len(set(pred_states))
    
    confusion = confusion_matrix(true_states, pred_states)
    
    # Pad matrix if needed
    if n_true != n_pred:
        max_dim = max(n_true, n_pred)
        padded_confusion = np.zeros((max_dim, max_dim))
        padded_confusion[:confusion.shape[0], :confusion.shape[1]] = confusion
        confusion = padded_confusion
    
    # Use Hungarian algorithm to find optimal assignment
    # We want to maximize overlap, so we use negative values
    row_ind, col_ind = linear_sum_assignment(-confusion)
    
    # Create mapping from predicted to true states
    state_mapping = {}
    for pred_idx, true_idx in zip(col_ind, row_ind):
        if pred_idx < n_pred and true_idx < n_true:
            state_mapping[pred_idx] = true_idx
    
    # Map predicted states
    aligned_pred_states = [state_mapping.get(s, s) for s in pred_states]
    
    return aligned_pred_states, state_mapping

def calculate_metrics(true_states, pred_states, aligned_pred_states):
    """Calculate various evaluation metrics"""
    
    metrics = {}
    
    # Basic metrics
    metrics['accuracy'] = accuracy_score(true_states, aligned_pred_states)
    
    # Clustering metrics (don't require alignment)
    metrics['adjusted_rand_score'] = adjusted_rand_score(true_states, pred_states)
    metrics['normalized_mutual_info'] = normalized_mutual_info_score(true_states, pred_states)
    
    # Precision, recall, F1 (macro average)
    precision, recall, f1, support = precision_recall_fscore_support(
        true_states, aligned_pred_states, average='macro', zero_division=0
    )
    metrics['precision_macro'] = precision
    metrics['recall_macro'] = recall
    metrics['f1_macro'] = f1
    
    # Per-class metrics
    precision_per_class, recall_per_class, f1_per_class, support_per_class = precision_recall_fscore_support(
        true_states, aligned_pred_states, average=None, zero_division=0
    )
    
    metrics['per_class'] = {}
    for i in range(len(precision_per_class)):
        metrics['per_class'][f'state_{i}'] = {
            'precision': float(precision_per_class[i]),
            'recall': float(recall_per_class[i]),
            'f1': float(f1_per_class[i]),
            'support': int(support_per_class[i])
        }
    
    # Confusion matrix
    cm = confusion_matrix(true_states, aligned_pred_states)
    metrics['confusion_matrix'] = cm.tolist()
    
    return metrics

def create_comparison_plots(df, results, aligned_pred_states, state_mapping, output_file):
    """Create comparison plots"""
    
    fig, axes = plt.subplots(4, 1, figsize=(15, 12))
    
    # Plot 1: Signal with true states
    colors_true = plt.cm.Set1(np.array(df['true_state']) / max(df['true_state']))
    axes[0].scatter(df['position'], df['signal'], c=colors_true, alpha=0.6, s=1)
    axes[0].set_ylabel('Signal')
    axes[0].set_title('Signal with True States')
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Signal with predicted states
    colors_pred = plt.cm.Set1(np.array(results['predicted_states']) / max(results['predicted_states']))
    axes[1].scatter(df['position'], df['signal'], c=colors_pred, alpha=0.6, s=1)
    axes[1].set_ylabel('Signal')
    axes[1].set_title('Signal with Predicted States (Original)')
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Signal with aligned predicted states
    colors_aligned = plt.cm.Set1(np.array(aligned_pred_states) / max(aligned_pred_states))
    axes[2].scatter(df['position'], df['signal'], c=colors_aligned, alpha=0.6, s=1)
    axes[2].set_ylabel('Signal')
    axes[2].set_title('Signal with Aligned Predicted States')
    axes[2].grid(True, alpha=0.3)
    
    # Plot 4: State comparison
    axes[3].plot(df['position'], df['true_state'], 'b-', alpha=0.7, linewidth=1, label='True States')
    axes[3].plot(df['position'], aligned_pred_states, 'r--', alpha=0.7, linewidth=1, label='Predicted States')
    axes[3].set_xlabel('Position')
    axes[3].set_ylabel('State')
    axes[3].set_title('State Comparison')
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def create_confusion_matrix_plot(confusion_matrix, output_file):
    """Create confusion matrix heatmap"""
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(confusion_matrix, annot=True, fmt='d', cmap='Blues', 
                xticklabels=[f'Pred {i}' for i in range(len(confusion_matrix[0]))],
                yticklabels=[f'True {i}' for i in range(len(confusion_matrix))])
    plt.title('Confusion Matrix (After State Alignment)')
    plt.ylabel('True State')
    plt.xlabel('Predicted State')
    plt.tight_layout()
    plt.savefig(output_file.replace('.png', '_confusion_matrix.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_posterior_probability_plot(results, df, output_file):
    """Create posterior probability visualization"""
    
    posterior_probs = np.array(results['posterior_probs'])
    n_states = posterior_probs.shape[1]
    
    fig, axes = plt.subplots(n_states + 1, 1, figsize=(15, 2 * (n_states + 1)))
    
    # Plot signal
    axes[0].plot(df['position'], df['signal'], 'k-', alpha=0.7, linewidth=0.5)
    axes[0].set_ylabel('Signal')
    axes[0].set_title('Original Signal')
    axes[0].grid(True, alpha=0.3)
    
    # Plot posterior probabilities for each state
    for i in range(n_states):
        axes[i+1].fill_between(df['position'], 0, posterior_probs[:, i], 
                              alpha=0.7, color=plt.cm.Set1(i))
        axes[i+1].set_ylabel(f'P(State {i})')
        axes[i+1].set_title(f'Posterior Probability - State {i}')
        axes[i+1].set_ylim(0, 1)
        axes[i+1].grid(True, alpha=0.3)
    
    axes[-1].set_xlabel('Position')
    plt.tight_layout()
    plt.savefig(output_file.replace('.png', '_posterior_probs.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generate_html_report(metrics, results, df, state_mapping, output_file):
    """Generate HTML evaluation report"""
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>HMM Segmentation Evaluation Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
            .metric {{ margin: 10px 0; }}
            .section {{ margin: 30px 0; }}
            .good {{ color: green; font-weight: bold; }}
            .warning {{ color: orange; font-weight: bold; }}
            .poor {{ color: red; font-weight: bold; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>HMM Segmentation Evaluation Report</h1>
            <p><strong>Dataset:</strong> {len(df)} genomic positions</p>
            <p><strong>Number of States:</strong> {results['model']['num_states']}</p>
            <p><strong>Processing Time:</strong> {results['processing_time_ms']} ms</p>
            <p><strong>Model Iterations:</strong> {results['model']['iterations']}</p>
            <p><strong>Final Log-likelihood:</strong> {results['model']['log_likelihood']:.6f}</p>
        </div>
        
        <div class="section">
            <h2>Overall Performance Metrics</h2>
    """
    
    # Add overall metrics
    def get_quality_class(value, thresholds):
        if value >= thresholds[0]:
            return "good"
        elif value >= thresholds[1]:
            return "warning"
        else:
            return "poor"
    
    accuracy_class = get_quality_class(metrics['accuracy'], [0.8, 0.6])
    ari_class = get_quality_class(metrics['adjusted_rand_score'], [0.7, 0.5])
    nmi_class = get_quality_class(metrics['normalized_mutual_info'], [0.7, 0.5])
    f1_class = get_quality_class(metrics['f1_macro'], [0.7, 0.5])
    
    html_content += f"""
            <div class="metric">
                <strong>Accuracy:</strong> <span class="{accuracy_class}">{metrics['accuracy']:.3f}</span>
            </div>
            <div class="metric">
                <strong>Adjusted Rand Score:</strong> <span class="{ari_class}">{metrics['adjusted_rand_score']:.3f}</span>
            </div>
            <div class="metric">
                <strong>Normalized Mutual Information:</strong> <span class="{nmi_class}">{metrics['normalized_mutual_info']:.3f}</span>
            </div>
            <div class="metric">
                <strong>F1 Score (Macro):</strong> <span class="{f1_class}">{metrics['f1_macro']:.3f}</span>
            </div>
        </div>
        
        <div class="section">
            <h2>Per-State Performance</h2>
            <table>
                <tr>
                    <th>State</th>
                    <th>Precision</th>
                    <th>Recall</th>
                    <th>F1 Score</th>
                    <th>Support</th>
                </tr>
    """
    
    # Add per-state metrics
    for state, state_metrics in metrics['per_class'].items():
        html_content += f"""
                <tr>
                    <td>{state}</td>
                    <td>{state_metrics['precision']:.3f}</td>
                    <td>{state_metrics['recall']:.3f}</td>
                    <td>{state_metrics['f1']:.3f}</td>
                    <td>{state_metrics['support']}</td>
                </tr>
        """
    
    html_content += """
            </table>
        </div>
        
        <div class="section">
            <h2>State Mapping</h2>
            <p>Predicted states were aligned with true states using the Hungarian algorithm:</p>
            <ul>
    """
    
    # Add state mapping
    for pred_state, true_state in state_mapping.items():
        html_content += f"<li>Predicted State {pred_state} → True State {true_state}</li>"
    
    html_content += """
            </ul>
        </div>
        
        <div class="section">
            <h2>Model Parameters</h2>
            <table>
                <tr>
                    <th>State</th>
                    <th>Emission Mean</th>
                    <th>Emission Variance</th>
                    <th>Initial Probability</th>
                </tr>
    """
    
    # Add model parameters
    for i, state_params in enumerate(results['model']['states']):
        html_content += f"""
                <tr>
                    <td>State {i}</td>
                    <td>{state_params['emission_mean']:.3f}</td>
                    <td>{state_params['emission_var']:.3f}</td>
                    <td>{state_params['initial_prob']:.3f}</td>
                </tr>
        """
    
    html_content += """
            </table>
        </div>
        
        <div class="section">
            <h2>Interpretation</h2>
            <h3>Metric Explanations:</h3>
            <ul>
                <li><strong>Accuracy:</strong> Percentage of correctly classified positions after state alignment</li>
                <li><strong>Adjusted Rand Score:</strong> Similarity between true and predicted clusterings (0-1, higher is better)</li>
                <li><strong>Normalized Mutual Information:</strong> Information shared between true and predicted states (0-1, higher is better)</li>
                <li><strong>F1 Score:</strong> Harmonic mean of precision and recall (0-1, higher is better)</li>
            </ul>
            
            <h3>Quality Guidelines:</h3>
            <ul>
                <li><span class="good">Good</span>: Accuracy ≥ 0.8, ARI/NMI ≥ 0.7, F1 ≥ 0.7</li>
                <li><span class="warning">Fair</span>: Accuracy ≥ 0.6, ARI/NMI ≥ 0.5, F1 ≥ 0.5</li>
                <li><span class="poor">Poor</span>: Below fair thresholds</li>
            </ul>
        </div>
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Evaluate HMM segmentation results')
    parser.add_argument('--results', required=True, help='Segmentation results JSON file')
    parser.add_argument('--truth', required=True, help='Ground truth CSV file')
    parser.add_argument('--output-report', default='evaluation_report.html',
                       help='Output HTML report file')
    parser.add_argument('--output-metrics', default='evaluation_metrics.json',
                       help='Output metrics JSON file')
    parser.add_argument('--output-plots', default='comparison_plots.png',
                       help='Output plots file')
    
    args = parser.parse_args()
    
    print("Loading segmentation results...")
    results = load_segmentation_results(args.results)
    
    print("Loading ground truth data...")
    df = load_truth_data(args.truth)
    
    print("Aligning predicted states with true states...")
    true_states = df['true_state'].tolist()
    pred_states = results['predicted_states']
    
    aligned_pred_states, state_mapping = align_states(true_states, pred_states)
    
    print("Calculating evaluation metrics...")
    metrics = calculate_metrics(true_states, pred_states, aligned_pred_states)
    
    print("Creating comparison plots...")
    create_comparison_plots(df, results, aligned_pred_states, state_mapping, args.output_plots)
    
    print("Creating confusion matrix plot...")
    create_confusion_matrix_plot(metrics['confusion_matrix'], args.output_plots)
    
    print("Creating posterior probability plots...")
    create_posterior_probability_plot(results, df, args.output_plots)
    
    print("Generating HTML report...")
    generate_html_report(metrics, results, df, state_mapping, args.output_report)
    
    print("Saving metrics...")
    # Add additional info to metrics
    metrics['state_mapping'] = state_mapping
    metrics['data_length'] = len(df)
    metrics['num_states'] = results['model']['num_states']
    metrics['processing_time_ms'] = results['processing_time_ms']
    metrics['model_iterations'] = results['model']['iterations']
    metrics['final_log_likelihood'] = results['model']['log_likelihood']
    
    with open(args.output_metrics, 'w') as f:
        json.dump(metrics, f, indent=2, default=str)
    
    print("\n=== Evaluation Summary ===")
    print(f"Accuracy: {metrics['accuracy']:.3f}")
    print(f"Adjusted Rand Score: {metrics['adjusted_rand_score']:.3f}")
    print(f"Normalized Mutual Information: {metrics['normalized_mutual_info']:.3f}")
    print(f"F1 Score (Macro): {metrics['f1_macro']:.3f}")
    print(f"Processing Time: {results['processing_time_ms']} ms")
    print(f"Model Iterations: {results['model']['iterations']}")
    
    print(f"\nState Mapping:")
    for pred_state, true_state in state_mapping.items():
        print(f"  Predicted State {pred_state} → True State {true_state}")
    
    print(f"\nOutputs:")
    print(f"  Report: {args.output_report}")
    print(f"  Metrics: {args.output_metrics}")
    print(f"  Plots: {args.output_plots}")

if __name__ == "__main__":
    main()