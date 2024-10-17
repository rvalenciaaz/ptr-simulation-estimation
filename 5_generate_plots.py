#!/usr/bin/env python3
# Usage: ./generate_plots.py <coverage_with_peaks.txt> <plot.png>

import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    data_file = sys.argv[1]
    plot_file = sys.argv[2]

    # Read data with peaks and troughs
    df = pd.read_csv(data_file, sep='\t')

    # Plot coverage
    plt.figure(figsize=(15, 5))
    plt.plot(df['pos'], df['coverage'], label='Coverage', color='gray')

    # Highlight peaks
    peaks = df[df['peak']]
    plt.scatter(peaks['pos'], peaks['coverage'], color='red', label='Peaks', s=10)

    # Highlight troughs
    troughs = df[df['trough']]
    plt.scatter(troughs['pos'], troughs['coverage'], color='blue', label='Troughs', s=10)

    plt.xlabel('Genomic Position')
    plt.ylabel('Coverage')
    plt.title('Coverage with Peaks and Troughs')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)

if __name__ == '__main__':
    main()
