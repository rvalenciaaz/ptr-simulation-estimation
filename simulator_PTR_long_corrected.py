import numpy as np
import matplotlib.pyplot as plt
import scipy.special

class Simulator:
    """Simulate long sequencing reads using the density of reads along a genome."""

    def __init__(self):
        # Define genome length
        self.genome_length = 4_000_000  # For example, E. coli genome length ~4 Mbp

        # Define bin edges and probabilities (uniform distribution)
        N = self.genome_length  # Probability for each position
        self.positions = np.arange(N)
        self.probs = np.ones(N) / N  # Uniform probabilities over positions

        self.min_log2_ptr = 0.05
        self.max_log2_ptr = 1.25

    def get_genome_length(self):
        """Return the number of bases in the genome."""
        return self.genome_length

    def simulate_reads(self, nreads, ptr=None, ori_pos=None):
        """Simulate long reads with specified ptr and origin position.

        Parameters
        ----------
            nreads : int
                Number of reads to simulate
            ptr : float
                PTR to simulate. If None, a PTR is chosen at random.
            ori_pos : float
                In the interval [0, 1]. Position of the replication origin.
                If None, a position is chosen at random.

        Returns
        -------
            reads : list of tuples
                List containing tuples of (start_position, end_position) for each read
            ptr : float
                The simulated PTR
            ori_pos : float
                The simulated replication origin as a fraction along the genome
            ter_pos : float
                The simulated replication terminus position as a fraction along the genome
            adj_probs : np.array
                Adjusted probabilities over the genome
        """
        if ptr is None:
            ptr = (
                np.random.random() * (self.max_log2_ptr - self.min_log2_ptr)
                + self.min_log2_ptr
            )
            ptr = 2 ** ptr
        else:
            assert np.log2(ptr) > self.min_log2_ptr, "ptr < 2^0.5 may be unreliable"

        if ori_pos is None:
            ori_pos = np.random.random()
        else:
            assert 0 <= ori_pos <= 1, "oriC coordinates must be in [0, 1]"

        ter_pos = (ori_pos + 0.5) % 1
        adj_probs = self.adjust_read_probs(ptr, ori_pos, ter_pos)

        # Simulate read start positions based on adjusted probabilities
        read_starts = np.random.choice(
            self.positions, size=nreads, p=adj_probs
        )

        # Simulate read lengths from a log-normal distribution typical for long reads
        read_lengths = self.simulate_read_lengths(nreads)

        # Ensure reads do not exceed genome boundaries (circular genome)
        reads = []
        for start, length in zip(read_starts, read_lengths):
            end = start + length
            if end < self.genome_length:
                reads.append((start, end))
            else:
                # Handle circular genome by wrapping around
                end = end % self.genome_length
                reads.append((start, end))

        return reads, ptr, ori_pos, ter_pos, adj_probs

    def simulate_read_lengths(self, nreads):
        """Simulate read lengths from a log-normal distribution.

        Parameters
        ----------
            nreads : int
                Number of reads to simulate

        Returns
        -------
            read_lengths : np.array
                Simulated read lengths
        """
        # Parameters for log-normal distribution (mean and sigma)
        mean_length = 10_000  # Mean read length of 10 kb
        sigma = 0.8  # Standard deviation in log space

        # Calculate mu parameter for log-normal distribution
        mu = np.log(mean_length) - 0.5 * sigma ** 2

        # Simulate read lengths
        read_lengths = np.random.lognormal(mean=mu, sigma=sigma, size=nreads)
        read_lengths = read_lengths.astype(int)

        # Cap maximum read length to genome length
        read_lengths = np.minimum(read_lengths, self.genome_length)

        return read_lengths

    def adjust_read_probs(self, ptr, ori_pos, ter_pos):
        """Scale position probabilities based on the PTR.

        Parameters
        ----------
            ptr : float
                The PTR to simulate
            ori_pos : float
                In [0,1]. The position of the replication origin
            ter_pos : float
                In [0,1]. The position of the replication terminus

        Returns
        -------
            adj_probs : np.array
                Adjusted probabilities over the genome positions
        """
        genome_length = self.genome_length
        positions = self.positions / genome_length  # Positions from 0 to 1

        # Compute replication profile based on PTR
        if ori_pos == ter_pos:
            raise ValueError("Origin and terminus positions cannot be the same.")

        alpha = np.log2(ptr) / (ori_pos - ter_pos)
        x1 = min(ori_pos, ter_pos)
        x2 = max(ori_pos, ter_pos)

        # Compute constants based on origin and terminus positions
        if ori_pos < ter_pos:
            c1 = np.log2(ptr)
            c2 = 0
        else:
            c1 = 0
            c2 = np.log2(ptr)

        # Adjust probabilities based on position
        m = positions
        log_probs = np.where(
            m <= x1,
            -alpha * (m - x1) + c1,
            np.where(
                m < x2,
                alpha * (m - x1) + c1,
                -alpha * (m - x2) + c2
            )
        )

        # Convert log probabilities to probabilities
        probs = 2 ** log_probs
        probs /= probs.sum()  # Normalize

        return probs

# Simulation and plotting
# -----------------------

# Instantiate the simulator
simulator = Simulator()

# Simulate reads
nreads = 2000  # Number of reads (adjust based on computational resources)
reads, ptr, ori_pos, ter_pos, adj_probs = simulator.simulate_reads(nreads)

# Compute coverage
genome_length = simulator.get_genome_length()
coverage = np.zeros(genome_length)

for start, end in reads:
    if start <= end:
        coverage[start:end] += 1
    else:
        # Handle circular genome wrapping
        coverage[start:] += 1
        coverage[:end] += 1

# Downsample coverage for plotting (e.g., bin the coverage)
bin_size = 10000  # Bin size of 10 kb
num_bins = genome_length // bin_size
binned_coverage = np.add.reduceat(coverage, np.arange(0, genome_length, bin_size)) / bin_size
bin_centers = np.arange(bin_size / 2, genome_length, bin_size)

# Expected coverage calculation
mean_read_length = 10_000  # As defined in simulate_read_lengths
total_bases = nreads * mean_read_length
expected_coverage_per_position = total_bases / genome_length
expected_coverage = expected_coverage_per_position * adj_probs * genome_length

# Binned expected coverage
expected_binned_coverage = np.add.reduceat(expected_coverage, np.arange(0, genome_length, bin_size)) / bin_size

# Plotting the simulated and expected coverage
plt.figure(figsize=(12, 6))
plt.plot(bin_centers, binned_coverage, label='Simulated coverage')
plt.plot(bin_centers, expected_binned_coverage, label='Expected coverage', linestyle='--')
plt.xlabel('Genome position (bp)')
plt.ylabel('Coverage (reads per base)')
plt.title('Simulated Long-Read Coverage Along the Genome')
plt.legend()
plt.tight_layout()
plt.savefig("simulated_vs_expected_long.png",dpi=800, bbox_inches="tight")
#plt.show()

# Plot read length distribution
read_lengths = [end - start if end >= start else (simulator.genome_length - start + end) for start, end in reads]
plt.figure(figsize=(8, 5))
plt.hist(read_lengths, bins=50, alpha=0.7)
plt.xlabel('Read Length (bp)')
plt.ylabel('Frequency')
plt.title('Simulated Long-Read Length Distribution')
plt.tight_layout()
plt.savefig("long_read_distribution.png",dpi=800, bbox_inches="tight")
#plt.show()

# Print simulated parameters
print(f"Simulated PTR: {ptr:.4f}")
print(f"Origin position (fraction of genome): {ori_pos:.4f}")
print(f"Terminus position (fraction of genome): {ter_pos:.4f}")