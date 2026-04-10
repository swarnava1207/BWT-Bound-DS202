import csv
import math
import time
from Bio import SeqIO

def compute_lz77_phrases(s):
    z, i, n = 0, 0, len(s)
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

def compute_sa_bwt_lcp(s):
    s = s + '$'
    n = len(s)
    sa = sorted(range(n), key=lambda i: s[i:])
    bwt = "".join(s[i - 1] if i > 0 else s[-1] for i in sa)
    
    lcp = [0] * n
    for i in range(1, n):
        idx1, idx2 = sa[i-1], sa[i]
        l = 0
        while idx1 + l < n and idx2 + l < n and s[idx1+l] == s[idx2+l]:
            l += 1
        lcp[i] = l
    return sa, bwt, lcp

def get_irreducible_lcps(bwt, lcp):
    irred_lcps = [lcp[0]]
    for i in range(1, len(bwt)):
        if bwt[i] != bwt[i-1]:
            irred_lcps.append(lcp[i])
    return irred_lcps

def generate_authors_string(length, alphabet_size=4):
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
    output_csv = "lemma31_lcp_powers_experiment.csv"
    string_length = 2000
    
    try:
        record = next(SeqIO.parse(fasta_file, "fasta"))
        dna_seq = str(record.seq).upper()[:string_length]
    except FileNotFoundError:
        print(f"Error: {fasta_file} not found.")
        return
        
    wc_seq = generate_authors_string(string_length, alphabet_size=4)
    
    z_dna = compute_lz77_phrases(dna_seq)
    _, bwt_dna, lcp_dna = compute_sa_bwt_lcp(dna_seq)
    irred_dna = get_irreducible_lcps(bwt_dna, lcp_dna)
    
    z_wc = compute_lz77_phrases(wc_seq)
    _, bwt_wc, lcp_wc = compute_sa_bwt_lcp(wc_seq)
    irred_wc = get_irreducible_lcps(bwt_wc, lcp_wc)
    
    bound_dna = 6 * z_dna * math.log2(string_length)
    bound_wc = 6 * z_wc * math.log2(string_length)
    
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["String_Type", "L_value", "Bucket_Range", "Irreducible_Count", "LZ77_z", "Theoretical_Max_Bound"])
        
        # Test l exclusively at powers of 2
        max_power = int(math.log2(string_length))
        l_values = [2**i for i in range(max_power)]
        
        for l in l_values:
            bucket_label = f"[{l}, {2*l})"
            
            count_dna = sum(1 for val in irred_dna if l <= val < 2*l)
            writer.writerow(["DNA", l, bucket_label, count_dna, z_dna, bound_dna])
            
            count_wc = sum(1 for val in irred_wc if l <= val < 2*l)
            writer.writerow(["Worst-Case", l, bucket_label, count_wc, z_wc, bound_wc])
                
    print(f"Finished! Data saved to {output_csv}")

if __name__ == "__main__":
    main()