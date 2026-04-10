import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

def main():
    csv_file = "lemma31_lcp_powers_experiment.csv"
    os.makedirs('res', exist_ok=True)
    
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: {csv_file} not found.")
        return
        
    df_dna = df[df['String_Type'] == 'DNA'].copy()
    df_wc = df[df['String_Type'] == 'Worst-Case'].copy()
    
    # Calculate the normalized ratio
    # Ratio = Count / (z * log_2(n))
    log2_n = math.log2(2000)
    df_dna['Normalized_Ratio'] = df_dna['Irreducible_Count'] / (df_dna['LZ77_z'] * log2_n)
    df_wc['Normalized_Ratio'] = df_wc['Irreducible_Count'] / (df_wc['LZ77_z'] * log2_n)
    
    # Setup for Bar Chart
    labels = df_wc['Bucket_Range'].tolist()
    x = np.arange(len(labels))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plotting the bars
    rects1 = ax.bar(x - width/2, df_wc['Normalized_Ratio'], width, label='Worst-Case Ratio', color='purple', alpha=0.8)
    rects2 = ax.bar(x + width/2, df_dna['Normalized_Ratio'], width, label='DNA Ratio', color='blue', alpha=0.8)
    
    # Add a line showing the theoretical max (6.0 based on the 6lz log n proof)
    ax.axhline(y=6.0, color='red', linestyle='--', linewidth=1.5, label='Proof Upper Bound Ceiling (6.0)')
    
    ax.set_xlabel('LCP Bucket Range $[\ell, 2\ell)$')
    ax.set_ylabel(r'Normalized Ratio: Count / $(z \log_2 n)$')
    ax.set_title('Distribution of Irreducible LCP Values Across Powers of 2 (Lemma 3.1)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45)
    ax.legend()
    
    plt.grid(axis='y', linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join('res', 'lemma31_binned_ratio.png'), dpi=300)
    plt.close()
    
    print("Plot saved to 'res/lemma31_binned_ratio.png'.")

if __name__ == "__main__":
    main()