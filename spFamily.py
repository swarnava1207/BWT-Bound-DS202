import csv
import time
import math

def compute_bwt_and_runs(s):
    """Computes the BWT of a string (appending '$') and counts the runs (r)."""
    s = s + '$'
    sa = sorted(range(len(s)), key=lambda i: s[i:])
    bwt = "".join(s[i - 1] if i > 0 else s[-1] for i in sa)
    
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
        for j in range(1, n - i + 1):
            if s.find(s[i:i+j], 0, i) != -1:
                max_len = j
            else:
                break
        i += max_len if max_len > 0 else 1
        z += 1
    return z

def generate_authors_string(length, alphabet_size):
    """
    Generates a substring of the authors' worst-case string T_{l,K}.
    Adapted for different alphabet sizes.
    """
    if alphabet_size == 2:
        base, sep = 2, ""
        K = 0 
    elif alphabet_size == 3:
        base, sep = 2, "2"
        K = max(1, int(math.log2(length) / 3))
    else: # alphabet_size == 4
        base, sep = 3, "3"
        K = max(1, int(math.log2(length) / 3))
        
    l = 2
    while (base**l) * (K + l) < length and l < 15:
        l += 1
        
    res = []
    for i in range(base**l):
        if base == 2:
            val = bin(i)[2:].zfill(l)
        else:
            n_val = i
            digits = []
            if n_val == 0:
                digits = ['0']
            while n_val:
                digits.append(str(n_val % base))
                n_val //= base
            val = "".join(digits[::-1]).zfill(l)
            
        res.append((sep * K) + val)
        if sum(len(b) for b in res) >= length:
            break
            
    return "".join(res)[:length]

def main():
    output_csv = "bwt_lz77_theoretical_bounds.csv"
    
    # Range parameters
    min_len = 10
    max_len = 4000
    alphabet_sizes = [2, 3, 4]
    batch_size = 50
    rows_buffer = []
    
    # Initialize CSV with only the raw, independent metrics
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "String_Length", "Alphabet_Size", "BWT_Runs_r", 
            "LZ77_Phrases_z", "BWT_String"
        ])

    total_combinations = (max_len - min_len + 1) * len(alphabet_sizes)
    print(f"Starting experiment: processing {total_combinations} theoretical strings...")
    start_time = time.time()
    processed = 0

    with open(output_csv, mode='a', newline='') as f:
        writer = csv.writer(f)
        
        for length in range(min_len, max_len + 1):
            for alpha_size in alphabet_sizes:
                
                # Generate string and compute raw metrics
                sub_seq = generate_authors_string(length, alpha_size)
                bwt_str, r = compute_bwt_and_runs(sub_seq)
                z = compute_lz77_phrases(sub_seq)
                
                rows_buffer.append([length, alpha_size, r, z, bwt_str])
                processed += 1
                
                # Save in batches
                if len(rows_buffer) >= batch_size:
                    writer.writerows(rows_buffer)
                    rows_buffer.clear()
                    
            if length % 200 == 0:
                print(f"Processed lengths up to {length}...")
                
        # Flush any remaining rows
        if rows_buffer:
            writer.writerows(rows_buffer)

    print(f"Finished! Processed {processed} strings in {time.time() - start_time:.2f} seconds.")
    print(f"Results saved to {output_csv}")

if __name__ == "__main__":
    main()