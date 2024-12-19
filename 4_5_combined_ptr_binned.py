#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # For running without a display
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv) < 5:
        print("Usage: ./combined_script.py <coverage.txt> <output.txt> <plot.png> <bin_size>")
        sys.exit(1)

    coverage_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_file = sys.argv[3]
    bin_size = int(sys.argv[4])

    # Read coverage data
    df = pd.read_csv(coverage_file, sep='\t', header=None, names=['chrom', 'pos', 'coverage'])

    # Bin coverage by the specified bin size
    # We assume a single chromosome or continuous coverage. If multiple chromosomes are present,
    # grouping by both chrom and bin_start would be advisable.
    df['bin_start'] = ( (df['pos'] - 1) // bin_size ) * bin_size + 1
    # Group by chromosome and bin start; compute mean coverage per bin
    binned = (df.groupby(['chrom', 'bin_start'], as_index=False)
                .agg({'coverage': 'mean'}))

    # Rename bin_start to pos for plotting convenience (the representative position of the bin)
    # Here we use the start of the bin as the representative position
    binned.rename(columns={'bin_start': 'pos'}, inplace=True)

    # Compute a rolling mean over the binned data
    # Note: The rolling window applies based on the number of rows, not genomic position distance.
    # If the bin size is large, consider adjusting the window accordingly.
    window_size = 1000  # can be adjusted as needed
    binned['rolling_mean'] = binned['coverage'].rolling(window=window_size, min_periods=1, center=False).mean()

    # Define thresholds for peaks and troughs
    mean_rolling = binned['rolling_mean'].mean()
    std_rolling = binned['rolling_mean'].std()
    peak_threshold = mean_rolling + std_rolling
    trough_threshold = mean_rolling - std_rolling

    # Identify peaks and troughs
    binned['peak'] = binned['rolling_mean'] > peak_threshold
    binned['trough'] = binned['rolling_mean'] < trough_threshold

    # Save the results (including peaks and troughs)
    # The output includes: chrom, pos (representative bin position), coverage, rolling_mean, peak, trough
    binned.to_csv(output_file, sep='\t', index=False)

    # Generate and save the plot
    plt.figure(figsize=(15, 5))
    plt.plot(binned['pos'], binned['coverage'], label='Binned Coverage', color='gray')

    peaks = binned[binned['peak']]
    troughs = binned[binned['trough']]

    if not peaks.empty:
        plt.scatter(peaks['pos'], peaks['coverage'], color='red', label='Peaks', s=10)
    if not troughs.empty:
        plt.scatter(troughs['pos'], troughs['coverage'], color='blue', label='Troughs', s=10)

    plt.xlabel('Genomic Position (binned)')
    plt.ylabel('Mean Coverage per Bin')
    plt.title('Coverage with Peaks and Troughs')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)

if __name__ == '__main__':
    main()
