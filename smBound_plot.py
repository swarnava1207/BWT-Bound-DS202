import os
import pandas as pd
import matplotlib.pyplot as plt

def main():
    csv_file = "smBound.csv"
    os.makedirs('res', exist_ok=True)
    
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: {csv_file} not found.")
        return
        
    df_dna = df[df['String_Type'] == 'DNA']
    df_wc = df[df['String_Type'] == 'Worst-Case']
    
    # Plot 1: Absolute Counts vs Theoretical Bounds
    plt.figure(figsize=(10, 6))
    
    # Worst Case
    plt.plot(df_wc['Window_m'], df_wc['Unique_Substrings_Sm'], color='purple', label='Worst-Case |Sm|', linewidth=2)
    plt.plot(df_wc['Window_m'], df_wc['Theoretical_Bound_mz'], color='purple', linestyle='--', label='Worst-Case Bound (m * z)')
    
    # DNA
    plt.plot(df_dna['Window_m'], df_dna['Unique_Substrings_Sm'], color='blue', label='DNA |Sm|', linewidth=2)
    plt.plot(df_dna['Window_m'], df_dna['Theoretical_Bound_mz'], color='blue', linestyle='--', label='DNA Bound (m * z)')
    
    plt.title(r'Lemma 3.1 Verification: $|\mathcal{S}_m|$ vs $m \cdot z$')
    plt.xlabel('Substring Window Length (m)')
    plt.ylabel('Count ($|\mathcal{S}_m|$)')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join('res', 'sm_bound.png'), dpi=300)
    plt.close()

    # Plot 2: Ratio of |Sm| / (m*z)
    plt.figure(figsize=(10, 6))
    plt.plot(df_wc['Window_m'], df_wc['Ratio'], color='purple', label='Worst-Case Ratio', linewidth=2)
    plt.plot(df_dna['Window_m'], df_dna['Ratio'], color='blue', label='DNA Ratio', linewidth=2)
    
    # Draw the hard limit at 1.0
    plt.axhline(y=1.0, color='red', linestyle='-', linewidth=1.5, label='Theoretical Max (1.0)')
    
    plt.title(r'Verification of Claim: $|\mathcal{S}_m| / (m \cdot z) <= 1$')
    plt.xlabel('Substring Window Length (m)')
    plt.ylabel(r'Ratio $|\mathcal{S}_m| / (m \cdot z)$')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join('res', 'sm_bound_ratio.png'), dpi=300)
    plt.close()
    
    print("Plots saved to 'res' directory.")

if __name__ == "__main__":
    main()