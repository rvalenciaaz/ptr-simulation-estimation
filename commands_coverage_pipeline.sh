

bash 1_map_reads_short.sh E4T11_R1_trimmed.fastq E4T11_R2_trimmed.fastq GCF_014131755.1_ASM1413175v1_genomic.fna output.sam 16
bash 2_sort_bam.sh output.sam output.bam
bash 3_compute_coverage.sh output.bam coverage_file.txt
