from airflow import DAG
from airflow.operators.bash import BashOperator
from datetime import datetime
import os

default_args = {
    'owner': 'airflow',
    'start_date': datetime(2023, 10, 16),
    'depends_on_past': False,
}

with DAG('peak_trough_pipeline', default_args=default_args, schedule_interval=None) as dag:

    # Define file paths (adjust paths as needed)
    base_dir = '/path/to/your/project'
    reads = os.path.join(base_dir, 'data', 'reads.fastq')
    reference = os.path.join(base_dir, 'data', 'reference.fasta')
    sam_file = os.path.join(base_dir, 'results', 'output.sam')
    sorted_bam = os.path.join(base_dir, 'results', 'sorted.bam')
    coverage_file = os.path.join(base_dir, 'results', 'coverage.txt')
    peak_trough_file = os.path.join(base_dir, 'results', 'coverage_with_peaks.txt')
    plot_file = os.path.join(base_dir, 'results', 'coverage_plot.png')

    # Ensure the results directory exists
    if not os.path.exists(os.path.join(base_dir, 'results')):
        os.makedirs(os.path.join(base_dir, 'results'))

    map_reads = BashOperator(
        task_id='map_reads',
        bash_command=f'./map_reads.sh {reads} {reference} {sam_file}',
        cwd=base_dir
    )

    sort_bam = BashOperator(
        task_id='sort_bam',
        bash_command=f'./sort_bam.sh {sam_file} {sorted_bam}',
        cwd=base_dir
    )

    compute_coverage = BashOperator(
        task_id='compute_coverage',
        bash_command=f'./compute_coverage.sh {sorted_bam} {coverage_file}',
        cwd=base_dir
    )

    peak_trough_calc = BashOperator(
        task_id='peak_trough_calculation',
        bash_command=f'./peak_trough_calculation.py {coverage_file} {peak_trough_file}',
        cwd=base_dir
    )

    generate_plots = BashOperator(
        task_id='generate_plots',
        bash_command=f'./generate_plots.py {peak_trough_file} {plot_file}',
        cwd=base_dir
    )

    # Define task dependencies
    map_reads >> sort_bam >> compute_coverage >> peak_trough_calc >> generate_plots
