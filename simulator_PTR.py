import numpy as np
import matplotlib.pyplot as plt
import scipy.special

class Simulator:
    """Simulate sequencing reads using the density of reads along a genome."""

    def __init__(self):
        # Define bin edges and probabilities (uniform distribution)
        N = 1000  # Number of bins
        self.bin_edges = np.linspace(0, 1, N + 1)
        self.probs = np.ones(N) / N  # Uniform probabilities over bins

        self.min_log2_ptr = 0.05
        self.max_log2_ptr = 1.25

    def get_genome_length(self):
        """Return the number of bases in the genome."""
        return 100 * self.probs.size

    def simulate_reads(self, nreads, ptr=None, ori_pos=None):
        """Simulate nreads along 100bp bins with specified ptr and
        origin position. If ptr or ori_pos are None, they are chosen
        randomly.

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
            read_positions : np.array
                Coordinates of simulated reads along the genome
            ptr : float
                The simulated PTR
            ori_pos : float
                The simulated replication origin as a fraction along the genome
            ter_pos : float
                The simulated replication terminus position as a fraction along the genome
            adj_probs : np.array
                Adjusted probabilities over bins
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
            assert ori_pos >= 0 and ori_pos <= 1, "oriC coordinates must be in [0, 1]"

        ter_pos = (ori_pos + 0.5) % 1
        adj_probs = self.adjust_read_probs(ptr, ori_pos, ter_pos)

        # Simulate reads based on adjusted probabilities
        binned_counts = np.random.multinomial(nreads, adj_probs)

        # Convert counts to read positions
        read_positions = []
        for i, c in enumerate(binned_counts):
            size = binned_counts.size
            m = 0.5 * (2 * i + 1) / size
            read_positions += [m for _ in range(c)]
        read_positions = self.get_genome_length() * np.array(read_positions)

        return read_positions, ptr, ori_pos, ter_pos, adj_probs

    def adjust_read_probs(self, ptr, ori_pos, ter_pos):
        """Scale bin probabilities based on the PTR.

        Parameters
        ----------
            ptr : float
                The PTR to simulate
            ori_pos : float
                In [0,1]. The position of the replication origin
            ter_pos : float
                In [0,1]. The position of the replication terminus
        """
        if ori_pos == ter_pos:
            raise ValueError("Origin and terminus positions cannot be the same.")

        alpha = np.log2(ptr) / (ori_pos - ter_pos)

        adj_probs = np.zeros(self.probs.size)
        for i, p in enumerate(self.probs):

            # Take the midpoint of each bin
            if i == self.bin_edges.size - 1:
                m = 0.5 * (self.bin_edges[i] + 1)
                length = 1 - self.bin_edges[i]
            else:
                m = 0.5 * (self.bin_edges[i] + self.bin_edges[i + 1])
                length = self.bin_edges[i + 1] - self.bin_edges[i]

            x1 = min(ori_pos, ter_pos)
            x2 = max(ori_pos, ter_pos)

            # Compute constants based on origin and terminus positions
            if ori_pos < ter_pos:
                c1 = np.log2(ptr)
                c2 = 0
            else:
                c1 = 0
                c2 = np.log2(ptr)

            if m <= x1:
                adj_probs[i] = -alpha * (m - x1) + c1
            elif x1 < m < x2:
                adj_probs[i] = alpha * (m - x1) + c1
            else:
                adj_probs[i] = -alpha * (m - x2) + c2
            adj_probs[i] += np.log2(length)

        # Normalize probabilities
        adj_probs -= np.log2(np.power(2, adj_probs).sum())
        adj_probs = adj_probs / np.log2(np.exp(1))

        # Reweight and normalize
        new_probs = np.zeros(self.probs.size)
        mask = self.probs != 0
        new_probs[mask] = (
            np.log(self.probs[mask]) + adj_probs[mask]
        )
        new_probs[mask] -= scipy.special.logsumexp(
            new_probs[mask]
        )
        new_probs[mask] = np.exp(new_probs[mask])

        return new_probs

# Simulation and plotting
# -----------------------

# Instantiate the simulator
simulator = Simulator()

# Simulate reads
nreads = 20000
read_positions, ptr, ori_pos, ter_pos, adj_probs = simulator.simulate_reads(nreads)

# Get genome length
L = simulator.get_genome_length()

# Get bin centers
N = simulator.probs.size
bin_edges = simulator.bin_edges * L  # Scale bin edges to genome length
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Get counts
counts, _ = np.histogram(read_positions, bins=bin_edges)

# Expected counts
expected_counts = nreads * adj_probs

# Plotting the simulated and expected read coverage
plt.figure(figsize=(10, 6))
plt.plot(bin_centers, counts, label='Simulated read counts')
plt.plot(bin_centers, expected_counts, label='Expected read counts', linestyle='--')
plt.xlabel('Genome position (bp)')
plt.ylabel('Read counts')
plt.title('Simulated Read Coverage Along the Genome')
plt.legend()
plt.tight_layout()
plt.savefig("simulated_vs_expected.png",dpi=800, bbox_inches="tight")
#plt.show()

# Print simulated parameters
print(f"Simulated PTR: {ptr:.4f}")
print(f"Origin position (fraction of genome): {ori_pos:.4f}")
print(f"Terminus position (fraction of genome): {ter_pos:.4f}")
