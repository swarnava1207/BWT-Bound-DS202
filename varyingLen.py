import random
import csv
import time
from Bio import SeqIO

def compute_bwt_and_runs(s):
    """Computes the BWT of a string (appending '$') and counts the runs (r)."""
    s = s + '$'
    
    # Fast Suffix Array construction for small to medium strings
    sa = sorted(range(len(s)), key=lambda i: s[i:])
    
    # Generate BWT from Suffix Array
    bwt = "".join(s[i - 1] if i > 0 else s[-1] for i in sa)
    
    # Count the number of runs (r) in the BWT
    runs = 1
    for i in range(1, len(bwt)):
        if bwt[i] != bwt[i-1]:
            runs += 1
            
    return bwt, runs

def compute_lz77_phrases(s):
    """Computes the number of phrases (z) in the LZ77 parsing of string s."""
    z = 0
    i = 0
    n = len(s)
    
    while i < n:
        max_len = 0
        # Find the longest prefix of s[i:] that occurs in s[:i]
        for j in range(1, n - i + 1):
            sub = s[i:i+j]
            if s.find(sub, 0, i) != -1:
                max_len = j
            else:
                break
                
        # If no previous factor starts here, phrase is a single character
        if max_len == 0:
            i += 1
        else:
            i += max_len
        z += 1
        
    return z

def main():
    fasta_file = "y_chromo.fasta"
    output_csv = "bwt_lz77_experiments.csv"
    
    print("Loading FASTA file...")
    try:
        record = next(SeqIO.parse(fasta_file, "fasta"))
        sequence = str(record.seq).upper()
    except FileNotFoundError:
        print(f"Error: Could not find {fasta_file}. Make sure it is in the same directory.")
        return

    seq_length = len(sequence)
    print(f"Loaded sequence of length: {seq_length}")

    # Experiment parameters
    min_len = 10
    max_len = 4000
    step = 15
    total_strings = 30000
    
    lengths = list(range(min_len, max_len + step, step))
    strings_per_length = total_strings // len(lengths)  # 10000 / 200 = 50
    
    batch_size = 50
    rows_buffer = []
    
    # Initialize CSV with headers
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["String_Length", "Start_Index", "BWT_Runs_r", "LZ77_Phrases_z", "BWT_String"])

    print("Starting experiment processing...")
    start_time = time.time()
    total_processed = 0

    # Process and append data
    with open(output_csv, mode='a', newline='') as f:
        writer = csv.writer(f)
        
        for length in lengths:
            for _ in range(strings_per_length):
                # Pick a random valid starting index
                start_idx = random.randint(0, seq_length - length)
                sub_seq = sequence[start_idx : start_idx + length]
                
                bwt_str, r = compute_bwt_and_runs(sub_seq)
                z = compute_lz77_phrases(sub_seq)
                
                rows_buffer.append([length, start_idx, r, z, bwt_str])
                total_processed += 1
                
                # Save every 20 rows
                if len(rows_buffer) >= batch_size:
                    writer.writerows(rows_buffer)
                    rows_buffer.clear()
                    
            print(f"Processed lengths up to {length}... ({total_processed}/{total_strings})")
            
        # Flush the remaining rows
        if rows_buffer:
            writer.writerows(rows_buffer)

    print(f"Finished! Processed {total_processed} strings in {time.time() - start_time:.2f} seconds.")
    print(f"Results safely saved to {output_csv}")

if __name__ == "__main__":
    main()