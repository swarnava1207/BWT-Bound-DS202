import random
import csv
import time
from Bio import SeqIO

def compute_bwt_runs(s):
    """Computes the BWT of a string (appending '$') and returns the run count."""
    s = s + '$'
    sa = sorted(range(len(s)), key=lambda i: s[i:])
    bwt = "".join(s[i - 1] if i > 0 else s[-1] for i in sa)
    
    runs = 1
    for i in range(1, len(bwt)):
        if bwt[i] != bwt[i-1]:
            runs += 1
    return runs


fasta_file = "y_chromo.fasta"
output_csv = "reverse_bwt_experiment.csv"

print("Loading FASTA file...")
try:
    record = next(SeqIO.parse(fasta_file, "fasta"))
    sequence = str(record.seq).upper()
except FileNotFoundError:
    print(f"Error: Could not find {fasta_file}.")

seq_length = len(sequence)

# Experiment Parameters
total_strings = 30000
min_len = 10
max_len = 4000
batch_size = 100
rows_buffer = []

with open(output_csv, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["String_Length", "Start_Index", "r_forward", "r_reverse"])

print(f"Starting experiment: processing {total_strings} strings...")
start_time = time.time()

with open(output_csv, mode='a', newline='') as f:
    writer = csv.writer(f)
    
    for i in range(total_strings):
        # Pick a random length and valid starting index
        length = random.randint(min_len, max_len)
        start_idx = random.randint(0, seq_length - length)
        
        # Extract forward string and create reverse string
        sub_seq_fwd = sequence[start_idx : start_idx + length]
        sub_seq_rev = sub_seq_fwd[::-1]
        
        # Compute runs for both
        r_fwd = compute_bwt_runs(sub_seq_fwd)
        r_rev = compute_bwt_runs(sub_seq_rev)
        
        rows_buffer.append([length, start_idx, r_fwd, r_rev])
        
        # Save in batches
        if len(rows_buffer) >= batch_size:
            writer.writerows(rows_buffer)
            rows_buffer.clear()
            
        if (i + 1) % 1000 == 0:
            print(f"Processed {i + 1}/{total_strings} strings...")
            
    # Flush remaining rows
    if rows_buffer:
        writer.writerows(rows_buffer)

print(f"Finished in {time.time() - start_time:.2f} seconds. Saved to {output_csv}")
