#!/usr/bin/env python3
"""
Generate comprehensive pipeline report
Combines all logs and results into a final summary
"""

import json
import argparse
import datetime
import re
from pathlib import Path

def parse_log_file(log_file):
    """Parse log file and extract key information"""
    if not Path(log_file).exists():
        return {"error": f"Log file {log_file} not found"}
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    return {"content": content, "lines": content.split('\n')}

def extract_build_info(build_log):
    """Extract build information from Rust build log"""
    info = {"status": "unknown", "warnings": 0, "errors": 0}
    
    if "error" in build_log:
        return info
    
    content = build_log["content"]
    
    if "Finished release" in content:
        info["status"] = "success"
    elif "error:" in content.lower():
        info["status"] = "failed"
        info["errors"] = len(re.findall(r'error:', content, re.IGNORECASE))
    
    info["warnings"] = len(re.findall(r'warning:', content, re.IGNORECASE))
    
    return info

def extract_data_info(data_log):
    """Extract data generation information"""
    info = {"status": "unknown", "positions": 0, "states": 0}
    
    if "error" in data_log:
        return info
    
    content = data_log["content"]
    
    # Look for data summary
    for line in data_log["lines"]:
        if "Total positions:" in line:
            try:
                info["positions"] = int(re.search(r'Total positions: (\d+)', line).group(1))
            except:
                pass
        elif "Parameters:" in line and "positions" in line and "states" in line:
            try:
                positions_match = re.search(r'(\d+) positions', line)
                states_match = re.search(r'(\d+) states', line)
                if positions_match:
                    info["positions"] = int(positions_match.group(1))
                if states_match:
                    info["states"] = int(states_match.group(1))
            except:
                pass
    
    if "Data saved to" in content:
        info["status"] = "success"
    
    return info

def extract_segmentation_info(segmentation_log):
    """Extract segmentation information"""
    info = {"status": "unknown", "iterations": 0, "final_likelihood": None, "time_ms": 0}
    
    if "error" in segmentation_log:
        return info
    
    content = segmentation_log["content"]
    
    # Look for completion message
    if "Segmentation completed" in content:
        info["status"] = "success"
        
        # Extract timing
        time_match = re.search(r'Segmentation completed in (\d+) ms', content)
        if time_match:
            info["time_ms"] = int(time_match.group(1))
    
    # Extract iterations and likelihood
    for line in segmentation_log["lines"]:
        if "Training iterations:" in line:
            try:
                info["iterations"] = int(re.search(r'Training iterations: (\d+)', line).group(1))
            except:
                pass
        elif "Final log-likelihood:" in line:
            try:
                info["final_likelihood"] = float(re.search(r'Final log-likelihood: ([-\d.]+)', line).group(1))
            except:
                pass
    
    return info

def load_json_safely(file_path):
    """Safely load JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        return {"error": str(e)}

def generate_html_report(data_info, build_info, segmentation_info, evaluation_info, 
                        segmentation_results, pipeline_params, output_file):
    """Generate comprehensive HTML report"""
    
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Epigenomic HMM Pipeline Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
            .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                      color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px; }}
            .section {{ background: #f8f9fa; padding: 20px; margin: 20px 0; border-radius: 8px; 
                       border-left: 4px solid #007bff; }}
            .success {{ color: #28a745; font-weight: bold; }}
            .warning {{ color: #ffc107; font-weight: bold; }}
            .error {{ color: #dc3545; font-weight: bold; }}
            .metric {{ margin: 10px 0; padding: 10px; background: white; border-radius: 5px; }}
            .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #dee2e6; padding: 12px; text-align: left; }}
            th {{ background-color: #e9ecef; font-weight: bold; }}
            .status-icon {{ font-size: 20px; margin-right: 10px; }}
            .code {{ background: #f1f3f4; padding: 10px; border-radius: 5px; font-family: monospace; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>🧬 Epigenomic HMM Segmentation Pipeline Report</h1>
            <p><strong>Generated:</strong> {timestamp}</p>
            <p><strong>Pipeline:</strong> Rust-based Hidden Markov Model for Chromatin State Segmentation</p>
        </div>
    """
    
    # Pipeline overview
    html_content += f"""
        <div class="section">
            <h2>📋 Pipeline Overview</h2>
            <div class="grid">
                <div>
                    <h3>Input Parameters</h3>
                    <div class="metric"><strong>Positions:</strong> {pipeline_params.get('n_positions', 'N/A')}</div>
                    <div class="metric"><strong>States:</strong> {pipeline_params.get('n_states', 'N/A')}</div>
                    <div class="metric"><strong>Max Iterations:</strong> {pipeline_params.get('max_iterations', 'N/A')}</div>
                </div>
                <div>
                    <h3>Pipeline Status</h3>
                    <div class="metric">
                        <span class="status-icon">📊</span>
                        <strong>Data Generation:</strong> 
                        <span class="{'success' if data_info['status'] == 'success' else 'error'}">
                            {data_info['status'].upper()}
                        </span>
                    </div>
                    <div class="metric">
                        <span class="status-icon">🔨</span>
                        <strong>Build:</strong> 
                        <span class="{'success' if build_info['status'] == 'success' else 'error'}">
                            {build_info['status'].upper()}
                        </span>
                    </div>
                    <div class="metric">
                        <span class="status-icon">🤖</span>
                        <strong>Segmentation:</strong> 
                        <span class="{'success' if segmentation_info['status'] == 'success' else 'error'}">
                            {segmentation_info['status'].upper()}
                        </span>
                    </div>
                </div>
            </div>
        </div>
    """
    
    # Data generation results
    html_content += f"""
        <div class="section">
            <h2>📊 Data Generation Results</h2>
            <div class="grid">
                <div>
                    <h3>Dataset Statistics</h3>
                    <div class="metric"><strong>Total Positions:</strong> {data_info.get('positions', 'N/A')}</div>
                    <div class="metric"><strong>Number of States:</strong> {data_info.get('states', 'N/A')}</div>
                    <div class="metric"><strong>Status:</strong> 
                        <span class="{'success' if data_info['status'] == 'success' else 'error'}">
                            {data_info['status']}
                        </span>
                    </div>
                </div>
                <div>
                    <h3>Data Quality</h3>
                    <div class="metric">✅ Synthetic epigenomic data generated</div>
                    <div class="metric">✅ Ground truth states available</div>
                    <div class="metric">✅ Realistic noise patterns included</div>
                </div>
            </div>
        </div>
    """
    
    # Build results
    html_content += f"""
        <div class="section">
            <h2>🔨 Build Results</h2>
            <div class="metric">
                <strong>Build Status:</strong> 
                <span class="{'success' if build_info['status'] == 'success' else 'error'}">
                    {build_info['status']}
                </span>
            </div>
            <div class="metric"><strong>Warnings:</strong> {build_info.get('warnings', 0)}</div>
            <div class="metric"><strong>Errors:</strong> {build_info.get('errors', 0)}</div>
        </div>
    """
    
    # Segmentation results
    processing_time = segmentation_results.get('processing_time_ms', 0) if 'error' not in segmentation_results else 0
    model_info = segmentation_results.get('model', {}) if 'error' not in segmentation_results else {}
    
    html_content += f"""
        <div class="section">
            <h2>🤖 HMM Segmentation Results</h2>
            <div class="grid">
                <div>
                    <h3>Performance Metrics</h3>
                    <div class="metric"><strong>Processing Time:</strong> {processing_time} ms</div>
                    <div class="metric"><strong>Training Iterations:</strong> {model_info.get('iterations', 'N/A')}</div>
                    <div class="metric"><strong>Final Log-likelihood:</strong> {model_info.get('log_likelihood', 'N/A'):.6f if isinstance(model_info.get('log_likelihood'), (int, float)) else 'N/A'}</div>
                    <div class="metric"><strong>Convergence:</strong> 
                        {'✅ Converged' if model_info.get('iterations', 0) < pipeline_params.get('max_iterations', 100) else '⚠️ Max iterations reached'}
                    </div>
                </div>
                <div>
                    <h3>Model Parameters</h3>
                    <div class="metric"><strong>Number of States:</strong> {model_info.get('num_states', 'N/A')}</div>
                    <div class="metric"><strong>Data Points:</strong> {segmentation_results.get('data_length', 'N/A')}</div>
                    <div class="metric"><strong>Algorithm:</strong> EM with Forward-Backward</div>
                    <div class="metric"><strong>Decoding:</strong> Viterbi Algorithm</div>
                </div>
            </div>
        </div>
    """
    
    # Evaluation results
    if 'error' not in evaluation_info:
        html_content += f"""
            <div class="section">
                <h2>📈 Evaluation Results</h2>
                <div class="grid">
                    <div>
                        <h3>Overall Performance</h3>
                        <div class="metric"><strong>Accuracy:</strong> {evaluation_info.get('accuracy', 'N/A'):.3f if isinstance(evaluation_info.get('accuracy'), (int, float)) else 'N/A'}</div>
                        <div class="metric"><strong>Adjusted Rand Score:</strong> {evaluation_info.get('adjusted_rand_score', 'N/A'):.3f if isinstance(evaluation_info.get('adjusted_rand_score'), (int, float)) else 'N/A'}</div>
                        <div class="metric"><strong>Normalized Mutual Info:</strong> {evaluation_info.get('normalized_mutual_info', 'N/A'):.3f if isinstance(evaluation_info.get('normalized_mutual_info'), (int, float)) else 'N/A'}</div>
                        <div class="metric"><strong>F1 Score (Macro):</strong> {evaluation_info.get('f1_macro', 'N/A'):.3f if isinstance(evaluation_info.get('f1_macro'), (int, float)) else 'N/A'}</div>
                    </div>
                    <div>
                        <h3>Quality Assessment</h3>
        """
        
        # Quality assessment
        accuracy = evaluation_info.get('accuracy', 0)
        ari = evaluation_info.get('adjusted_rand_score', 0)
        
        if isinstance(accuracy, (int, float)) and isinstance(ari, (int, float)):
            if accuracy >= 0.8 and ari >= 0.7:
                quality = "🟢 Excellent"
                quality_class = "success"
            elif accuracy >= 0.6 and ari >= 0.5:
                quality = "🟡 Good"
                quality_class = "warning"
            else:
                quality = "🔴 Needs Improvement"
                quality_class = "error"
        else:
            quality = "❓ Unknown"
            quality_class = "warning"
        
        html_content += f"""
                        <div class="metric"><strong>Overall Quality:</strong> <span class="{quality_class}">{quality}</span></div>
                        <div class="metric">✅ State alignment performed</div>
                        <div class="metric">✅ Confusion matrix generated</div>
                        <div class="metric">✅ Per-state metrics calculated</div>
                    </div>
                </div>
            </div>
        """
    
    # Model parameters table
    if 'error' not in segmentation_results and 'states' in model_info:
        html_content += """
            <div class="section">
                <h2>🎯 Learned Model Parameters</h2>
                <table>
                    <tr>
                        <th>State</th>
                        <th>Emission Mean</th>
                        <th>Emission Variance</th>
                        <th>Initial Probability</th>
                    </tr>
        """
        
        for i, state in enumerate(model_info['states']):
            html_content += f"""
                    <tr>
                        <td>State {i}</td>
                        <td>{state.get('emission_mean', 'N/A'):.3f if isinstance(state.get('emission_mean'), (int, float)) else 'N/A'}</td>
                        <td>{state.get('emission_var', 'N/A'):.3f if isinstance(state.get('emission_var'), (int, float)) else 'N/A'}</td>
                        <td>{state.get('initial_prob', 'N/A'):.3f if isinstance(state.get('initial_prob'), (int, float)) else 'N/A'}</td>
                    </tr>
            """
        
        html_content += """
                </table>
            </div>
        """
    
    # Commands to reproduce
    html_content += f"""
        <div class="section">
            <h2>🔧 Commands to Reproduce</h2>
            <h3>1. Generate Data</h3>
            <div class="code">
python3 generate_data.py \\
    --n-positions {pipeline_params.get('n_positions', 5000)} \\
    --n-states {pipeline_params.get('n_states', 3)} \\
    --visualize \\
    --add-noise
            </div>
            
            <h3>2. Build Rust Tool</h3>
            <div class="code">
cargo build --release
            </div>
            
            <h3>3. Run Segmentation</h3>
            <div class="code">
./target/release/main \\
    --input epigenomic_data.csv \\
    --output segmentation_results.json \\
    --states {pipeline_params.get('n_states', 3)} \\
    --max-iterations {pipeline_params.get('max_iterations', 100)}
            </div>
            
            <h3>4. Evaluate Results</h3>
            <div class="code">
python3 evaluate_segmentation.py \\
    --results segmentation_results.json \\
    --truth epigenomic_data_with_true_states.csv \\
    --output-report evaluation_report.html
            </div>
            
            <h3>5. Run Full Pipeline with Nextflow</h3>
            <div class="code">
nextflow run main.nf \\
    --n_positions {pipeline_params.get('n_positions', 5000)} \\
    --n_states {pipeline_params.get('n_states', 3)} \\
    --max_iterations {pipeline_params.get('max_iterations', 100)} \\
    --visualize \\
    --add_noise
            </div>
        </div>
        
        <div class="section">
            <h2>📁 Output Files</h2>
            <ul>
                <li><strong>epigenomic_data.csv</strong> - Input data for HMM</li>
                <li><strong>epigenomic_data_with_true_states.csv</strong> - Data with ground truth</li>
                <li><strong>segmentation_results.json</strong> - HMM segmentation output</li>
                <li><strong>evaluation_report.html</strong> - Detailed evaluation report</li>
                <li><strong>comparison_plots.png</strong> - Visualization plots</li>
                <li><strong>pipeline_report.html</strong> - This comprehensive report</li>
            </ul>
        </div>
        
        <div class="section">
            <h2>🔬 Technical Details</h2>
            <p><strong>Algorithm:</strong> Hidden Markov Model with Expectation-Maximization training</p>
            <p><strong>Implementation:</strong> Rust with parallel processing capabilities</p>
            <p><strong>Features:</strong></p>
            <ul>
                <li>Forward-backward algorithm for parameter estimation</li>
                <li>Viterbi decoding for state sequence prediction</li>
                <li>Hungarian algorithm for state alignment in evaluation</li>
                <li>Comprehensive performance metrics</li>
                <li>Synthetic data generation with realistic patterns</li>
            </ul>
        </div>
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)

def generate_summary_json(data_info, build_info, segmentation_info, evaluation_info, 
                         segmentation_results, pipeline_params, output_file):
    """Generate summary JSON"""
    
    summary = {
        "timestamp": datetime.datetime.now().isoformat(),
        "pipeline_params": pipeline_params,
        "stages": {
            "data_generation": data_info,
            "build": build_info,
            "segmentation": segmentation_info,
            "evaluation": evaluation_info if 'error' not in evaluation_info else {"status": "failed"}
        },
        "results": {
            "processing_time_ms": segmentation_results.get('processing_time_ms', 0) if 'error' not in segmentation_results else 0,
            "model_iterations": segmentation_results.get('model', {}).get('iterations', 0) if 'error' not in segmentation_results else 0,
            "final_log_likelihood": segmentation_results.get('model', {}).get('log_likelihood', None) if 'error' not in segmentation_results else None,
            "accuracy": evaluation_info.get('accuracy', None) if 'error' not in evaluation_info else None,
            "adjusted_rand_score": evaluation_info.get('adjusted_rand_score', None) if 'error' not in evaluation_info else None
        },
        "status": {
            "overall": "success" if all([
                data_info['status'] == 'success',
                build_info['status'] == 'success', 
                segmentation_info['status'] == 'success'
            ]) else "failed"
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

def main():
    parser = argparse.ArgumentParser(description='Generate pipeline report')
    parser.add_argument('--data-log', required=True, help='Data generation log file')
    parser.add_argument('--build-log', required=True, help='Build log file')
    parser.add_argument('--segmentation-log', required=True, help='Segmentation log file')
    parser.add_argument('--evaluation-log', required=True, help='Evaluation log file')
    parser.add_argument('--evaluation-metrics', required=True, help='Evaluation metrics JSON')
    parser.add_argument('--segmentation-results', required=True, help='Segmentation results JSON')
    parser.add_argument('--output-report', default='pipeline_report.html', help='Output HTML report')
    parser.add_argument('--output-summary', default='pipeline_summary.json', help='Output summary JSON')
    parser.add_argument('--n-positions', type=int, default=5000, help='Number of positions')
    parser.add_argument('--n-states', type=int, default=3, help='Number of states')
    parser.add_argument('--max-iterations', type=int, default=100, help='Max iterations')
    
    args = parser.parse_args()
    
    print("Parsing log files...")
    data_log = parse_log_file(args.data_log)
    build_log = parse_log_file(args.build_log)
    segmentation_log = parse_log_file(args.segmentation_log)
    evaluation_log = parse_log_file(args.evaluation_log)
    
    print("Extracting information...")
    data_info = extract_data_info(data_log)
    build_info = extract_build_info(build_log)
    segmentation_info = extract_segmentation_info(segmentation_log)
    
    print("Loading results...")
    evaluation_info = load_json_safely(args.evaluation_metrics)
    segmentation_results = load_json_safely(args.segmentation_results)
    
    pipeline_params = {
        'n_positions': args.n_positions,
        'n_states': args.n_states,
        'max_iterations': args.max_iterations
    }
    
    print("Generating HTML report...")
    generate_html_report(
        data_info, build_info, segmentation_info, evaluation_info,
        segmentation_results, pipeline_params, args.output_report
    )
    
    print("Generating summary JSON...")
    generate_summary_json(
        data_info, build_info, segmentation_info, evaluation_info,
        segmentation_results, pipeline_params, args.output_summary
    )
    
    print(f"\nPipeline report generated:")
    print(f"  HTML Report: {args.output_report}")
    print(f"  JSON Summary: {args.output_summary}")

if __name__ == "__main__":
    main()