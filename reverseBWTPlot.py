import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    csv_file = "reverse_bwt_experiment.csv"
    os.makedirs('res', exist_ok=True)
    
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: {csv_file} not found.")
        return
        
    # Calculate the bounding metric: r_reverse / (r_forward * log2(n)^2)
    df['log2_n_sq'] = np.log2(df['String_Length']) ** 2
    df['normalized_ratio'] = df['r_reverse'] / (df['r_forward'] * df['log2_n_sq'])
    
    # Plot 1: Direct Correlation between r and r_reverse
    plt.figure(figsize=(10, 6))
    plt.scatter(df['r_forward'], df['r_reverse'], alpha=0.3, color='teal', s=10)
    
    # Add a y=x reference line
    max_val = max(df['r_forward'].max(), df['r_reverse'].max())
    plt.plot([0, max_val], [0, max_val], color='red', linestyle='--', label='y = x')
    
    plt.title(r'BWT Runs: Forward ($r$) vs Reverse ($\bar{r}$)')
    plt.xlabel(r'Forward Runs ($r$)')
    plt.ylabel(r'Reverse Runs ($\bar{r}$)')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join('res', 'corollary311_r_vs_r_rev.png'), dpi=300)
    plt.close()

    # Plot 2: Verification of the Upper Bound
    plt.figure(figsize=(10, 6))
    plt.scatter(df['String_Length'], df['normalized_ratio'], alpha=0.2, color='darkgreen', s=10)
    
    plt.title(r'Verification of Corollary 3.11: $\bar{r} / (r \log_2^2 n) <= \mathcal{O}(1)$')
    plt.xlabel('String Length ($n$)')
    plt.ylabel(r'Normalized Ratio: $\bar{r} / (r \log_2^2 n)$')
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join('res', 'corollary311_bound.png'), dpi=300)
    plt.close()
    
    print("Plots successfully saved to 'res/corollary311_r_vs_r_rev.png' and 'res/corollary311_bound.png'.")

if __name__ == "__main__":
    main()