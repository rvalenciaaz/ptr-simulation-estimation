#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_genome_file> <chunk_size> <output_folder>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    chunk_size = int(sys.argv[2])
    output_folder = sys.argv[3]

    # Create the output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    # Read all sequences from the input file and concatenate them
    sequences = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequences.append(str(record.seq))
    full_sequence = "".join(sequences)
    seq_length = len(full_sequence)

    # Calculate how many chunks we need
    num_chunks = (seq_length + chunk_size - 1) // chunk_size

    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, seq_length)
        chunk_seq = full_sequence[start:end]

        # Create a Biopython SeqRecord for the chunk
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        chunk_record = SeqRecord(Seq(chunk_seq),
                                 id=f"chunk_{i+1}",
                                 description="")

        # Write the chunk to a separate FASTA file
        chunk_filename = os.path.join(output_folder, f"chunk_{i+1}.fa")
        SeqIO.write(chunk_record, chunk_filename, "fasta")

    print(f"Successfully split {seq_length} bases into {num_chunks} chunks.")

if __name__ == "__main__":
    main()
