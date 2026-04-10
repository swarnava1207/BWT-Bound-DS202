import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def compute_metrics(df):
    """Computes the derived metrics for plotting."""
    # Filter out any rows where z is 0 to avoid division by zero
    df = df[df['LZ77_Phrases_z'] > 0].copy()
    
    n = df['String_Length']
    r = df['BWT_Runs_r']
    z = df['LZ77_Phrases_z']
    
    df['r_over_z'] = r / z
    # Compute r / (z * log_2(n)^2)
    df['z_log2_n_sq'] = z * (np.log2(n) ** 2)
    df['bound_ratio'] = r / df['z_log2_n_sq']
    
    return df

def plot_and_save(datasets, x_col, y_col, xlabel, ylabel, title, filename):
    """Plots multiple datasets on the same axes and saves the figure."""
    plt.figure(figsize=(10, 6))
    
    for label, df, color, marker, alpha in datasets:
        if x_col in df.columns and y_col in df.columns:
            plt.scatter(df[x_col], df[y_col], label=label, color=color, marker=marker, alpha=alpha, s=15)
            
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join('res', filename), dpi=300)
    plt.close()

def main():
    # 1. Create results directory
    os.makedirs('res', exist_ok=True)
    print("Created 'res' directory.")

    # 2. Load the datasets
    try:
        df_var = pd.read_csv('varyingLen.csv')
        df_same = pd.read_csv('sameLen.csv')
        df_sp = pd.read_csv('spFamily.csv')
    except FileNotFoundError as e:
        print(f"Error loading files: {e}")
        return

    # 3. Compute metrics for all dataframes
    df_var = compute_metrics(df_var)
    df_same = compute_metrics(df_same)
    df_sp = compute_metrics(df_sp)

    # 4. Split the theoretical family by alphabet size for clearer coloring
    df_sp_2 = df_sp[df_sp['Alphabet_Size'] == 2]
    df_sp_3 = df_sp[df_sp['Alphabet_Size'] == 3]
    df_sp_4 = df_sp[df_sp['Alphabet_Size'] == 4]

    # Define the datasets to plot together: (Label, DataFrame, Color, Marker, Alpha)
    plot_groups = [
        ('DNA (Varying Len)', df_var, 'blue', 'o', 0.3),
        ('DNA (Len=1000)', df_same, 'green', 'x', 0.1),
        ('Worst-Case (Alpha=2)', df_sp_2, 'red', '^', 0.6),
        ('Worst-Case (Alpha=3)', df_sp_3, 'orange', 's', 0.6),
        ('Worst-Case (Alpha=4)', df_sp_4, 'purple', 'D', 0.6)
    ]

    print("Generating plots...")

    # Plot 1: r vs n
    plot_and_save(plot_groups, 'String_Length', 'BWT_Runs_r', 
                  'String Length (n)', 'BWT Runs (r)', 
                  'Growth of BWT Runs vs String Length', 'r_vs_n.png')

    # Plot 2: z vs n
    plot_and_save(plot_groups, 'String_Length', 'LZ77_Phrases_z', 
                  'String Length (n)', 'LZ77 Phrases (z)', 
                  'Growth of LZ77 Phrases vs String Length', 'z_vs_n.png')

    # Plot 3: r vs z
    plot_and_save(plot_groups, 'LZ77_Phrases_z', 'BWT_Runs_r', 
                  'LZ77 Phrases (z)', 'BWT Runs (r)', 
                  'BWT Runs (r) vs LZ77 Phrases (z)', 'r_vs_z.png')

    # Plot 4: r/z vs n
    plot_and_save(plot_groups, 'String_Length', 'r_over_z', 
                  'String Length (n)', 'Ratio (r/z)', 
                  'Ratio of r/z vs String Length', 'r_over_z_vs_n.png')

    # Plot 5: r / (z * log^2 n) vs n (The theoretical bound check)
    plot_and_save(plot_groups, 'String_Length', 'bound_ratio', 
                  'String Length (n)', r'$r / (z \log_2^2 n)$', 
                  'Verification of the BWT Conjecture Upper Bound', 'bound_verification.png')
    
    # Plot 6: r / (z * log^2 n) vs n (but just show alpha = 4 and DNA (Varying Len))
    plot_and_save([
        ('DNA (Varying Len)', df_var, 'blue', 'o', 0.3),
        ('Worst-Case (Alpha=4)', df_sp_4, 'purple', 'D', 0.6)
    ], 'String_Length', 'bound_ratio', 
       'String Length (n)', r'$r / (z \log_2^2 n)$', 
       'Verification of the BWT Conjecture Upper Bound (Alpha=4 vs DNA)', 'bound_verification_alpha4.png')

    print("All plots successfully saved to the 'res' directory.")

if __name__ == "__main__":
    main()