#!/usr/bin/env python3
# Usage: ./peak_trough_calculation.py <coverage.txt> <output.txt>

import sys
import pandas as pd

def main():
    coverage_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read coverage data
    df = pd.read_csv(coverage_file, sep='\t', header=None, names=['chrom', 'pos', 'coverage'])

    # Calculate rolling mean to smooth coverage data
    df['rolling_mean'] = df['coverage'].rolling(window=1000, min_periods=1).mean()

    # Define thresholds for peaks and troughs
    peak_threshold = df['rolling_mean'].mean() + df['rolling_mean'].std()
    trough_threshold = df['rolling_mean'].mean() - df['rolling_mean'].std()

    # Identify peaks and troughs
    df['peak'] = df['rolling_mean'] > peak_threshold
    df['trough'] = df['rolling_mean'] < trough_threshold

    # Save the results
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()
