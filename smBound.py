import csv
import math
import time
from Bio import SeqIO

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

def count_unique_cyclic_substrings(s, m):
    """
    Counts the number of unique length-m substrings in T^infinity.
    We simulate T^infinity by appending the first m-1 characters to the end.
    """
    if m == 0: return 0
    # Create the cyclic sequence just long enough to capture all length-m windows
    cyclic_s = s + s[:m-1]
    seen = set()
    for i in range(len(s)):
        seen.add(cyclic_s[i:i+m])
    return len(seen)

def generate_authors_string(length, alphabet_size=4):
    """Generates the Alpha=4 worst-case string from the paper."""
    base, sep = 3, "3"
    K = max(1, int(math.log2(length) / 3))
    l = 2
    while (base**l) * (K + l) < length and l < 15:
        l += 1
        
    res = []
    for i in range(base**l):
        n_val = i
        digits = []
        if n_val == 0: digits = ['0']
        while n_val:
            digits.append(str(n_val % base))
            n_val //= base
        val = "".join(digits[::-1]).zfill(l)
        res.append((sep * K) + val)
        if sum(len(b) for b in res) >= length:
            break
    return "".join(res)[:length]

def main():
    fasta_file = "y_chromo.fasta"
    output_csv = "smBound.csv"
    string_length = 4000
    max_window_m = 2000 # Test substring lengths from 1 to 2000
    
    # 1. Get DNA String
    try:
        record = next(SeqIO.parse(fasta_file, "fasta"))
        dna_seq = str(record.seq).upper()[:string_length]
    except FileNotFoundError:
        print(f"Error: {fasta_file} not found.")
        return
        
    # 2. Get Worst-Case String
    worst_case_seq = generate_authors_string(string_length, alphabet_size=4)
    
    # 3. Compute z for both
    z_dna = compute_lz77_phrases(dna_seq)
    z_wc = compute_lz77_phrases(worst_case_seq)
    
    print(f"DNA z = {z_dna} | Worst-Case z = {z_wc}")
    
    # 4. Sweep window lengths and record
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["String_Type", "Window_m", "Unique_Substrings_Sm", "LZ77_z", "Theoretical_Bound_mz", "Ratio"])
        
        start_time = time.time()
        for m in range(1, max_window_m + 1):
            # DNA
            sm_dna = count_unique_cyclic_substrings(dna_seq, m)
            bound_dna = m * z_dna
            writer.writerow(["DNA", m, sm_dna, z_dna, bound_dna, sm_dna / bound_dna if bound_dna > 0 else 0])
            
            # Worst-Case
            sm_wc = count_unique_cyclic_substrings(worst_case_seq, m)
            bound_wc = m * z_wc
            writer.writerow(["Worst-Case", m, sm_wc, z_wc, bound_wc, sm_wc / bound_wc if bound_wc > 0 else 0])
            
            if m % 100 == 0:
                print(f"Processed window size m={m}...")
                
    print(f"Finished in {time.time() - start_time:.2f} seconds. Saved to {output_csv}")

if __name__ == "__main__":
    main()