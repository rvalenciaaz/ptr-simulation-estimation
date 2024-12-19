#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # For running without a display
import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Bin coverage data, identify peaks and troughs, calculate PTR (peak-to-trough ratio), "
                    "and plot log2(PTR) along with rolling mean coverage."
    )
    parser.add_argument("coverage_file", type=str,
                        help="Input coverage file with columns: chrom, pos, coverage.")
    parser.add_argument("output_file", type=str,
                        help="Output TSV file with binning, PTR, peak, and trough annotations.")
    parser.add_argument("plot_file", type=str,
                        help="Output plot file (e.g., PNG) showing log2(PTR), peaks, troughs, and rolling mean.")
    parser.add_argument("bin_size", type=int,
                        help="Bin size (in bp) for grouping coverage data, e.g. 100.")
    parser.add_argument("window_size", type=int,
                        help="Window size for rolling mean calculation, e.g. 1000.")

    args = parser.parse_args()

    coverage_file = args.coverage_file
    output_file = args.output_file
    plot_file = args.plot_file
    bin_size = args.bin_size
    window_size = args.window_size

    # Read coverage data
    df = pd.read_csv(coverage_file, sep='\t', header=None, names=['chrom', 'pos', 'coverage'])

    # Bin coverage by the specified bin size
    df['bin_start'] = ((df['pos'] - 1) // bin_size) * bin_size + 1
    binned = (df.groupby(['chrom', 'bin_start'], as_index=False)
                .agg({'coverage': 'mean'}))
    binned.rename(columns={'bin_start': 'pos'}, inplace=True)

    # Compute a rolling mean over the binned data
    binned['rolling_mean'] = binned['coverage'].rolling(window=window_size, min_periods=1, center=False).mean()

    # Define thresholds for peaks and troughs based on rolling mean
    mean_rolling = binned['rolling_mean'].mean()
    std_rolling = binned['rolling_mean'].std()
    peak_threshold = mean_rolling + std_rolling
    trough_threshold = mean_rolling - std_rolling

    # Identify peaks and troughs
    binned['peak'] = binned['rolling_mean'] > peak_threshold
    binned['trough'] = binned['rolling_mean'] < trough_threshold

    # Compute PTR and log2(PTR)
    binned['PTR'] = binned['coverage'] / binned['rolling_mean']
    binned['log2_PTR'] = np.log2(binned['PTR'])

    # Save the results
    binned.to_csv(output_file, sep='\t', index=False)

    # Plotting
    fig, ax1 = plt.subplots(figsize=(15,5))

    # Plot log2(PTR) on the primary axis
    ax1.plot(binned['pos'], binned['log2_PTR'], label='log2(PTR)', color='gray')
    ax1.set_xlabel('Genomic Position (binned)')
    ax1.set_ylabel('log2(PTR)', color='gray')
    ax1.tick_params(axis='y', labelcolor='gray')

    # Scatter peaks and troughs on the log2(PTR) axis
    peaks = binned[binned['peak']]
    troughs = binned[binned['trough']]

    if not peaks.empty:
        ax1.scatter(peaks['pos'], peaks['log2_PTR'], color='red', label='Peaks', s=10)
    if not troughs.empty:
        ax1.scatter(troughs['pos'], troughs['log2_PTR'], color='blue', label='Troughs', s=10)

    # Create a secondary y-axis for rolling mean coverage
    ax2 = ax1.twinx()
    ax2.plot(binned['pos'], binned['rolling_mean'], label='Rolling Mean Coverage', color='green', alpha=0.6)
    ax2.set_ylabel('Rolling Mean Coverage', color='green')
    ax2.tick_params(axis='y', labelcolor='green')

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    plt.title('Coverage as log2(PTR) with Peaks and Troughs + Rolling Mean Coverage')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)

if __name__ == '__main__':
    main()
