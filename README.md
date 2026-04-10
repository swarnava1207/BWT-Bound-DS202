# BWT Runs and LZ77 Parsing Bounds Validation

This repository contains Python scripts designed to empirically test and visualize the theoretical relationships and bounds between the number of runs in the Burrows-Wheeler Transform ($r$) and the number of phrases in LZ77 parsing ($z$). 

The experiments compare real-world DNA data (extracted from `y_chromo.fasta`) against a theoretical "worst-case" string family described in the literature.

## Repository Structure

The repository is divided into scripts that generate data (CSVs) and scripts that generate visualizations (PNGs in the `res/` folder).

### 1. Main Relationship Experiments ($r$ vs $z$)
These scripts explore the general growth and ratio of $r$ and $z$ as string length ($n$) increases.

* **`varyingLen.py`**: Samples random sub-sequences of varying lengths (from 10 to 4000) from the FASTA file. Computes $r$ and $z$.
    * *Outputs:* `varyingLen.csv`
* **`samLen.py`**: Samples 20,000 random strings of a *fixed* length (2000) from the FASTA file to see the variance of $r$ and $z$ at a constant $n$.
    * *Outputs:* `sameLen.csv`
* **`spFamily.py`**: Generates the "worst-case" string family from the theoretical paper for alphabet sizes 2, 3, and 4. Computes $r$ and $z$ for these strings across various lengths.
    * *Outputs:* `spFamily.csv`
* **`plots.py`**: The master plotting script for the main experiments. It reads the three CSVs above and generates multiple comparisons.

### 2. Lemma & Corollary Verifications
These scripts test specific mathematical lemmas and corollaries related to the $r$ bounded by $z \log^2 n$ proof.

* **`smBound.py`**: Verifies Lemma 3.1 regarding the number of unique length-$m$ cyclic substrings ($|S_m|$). It tests the claim that $|S_m| \le m \cdot z$.
    * *Outputs:* `smBound.csv`
* **`smBound_plot.py`**: Reads `smBound.csv` to visualize the substring bounds.
* **`lem3_1_lcp.py`**: Another facet of Lemma 3.1, examining the distribution of irreducible Longest Common Prefix (LCP) values binned into powers of 2.
    * *Outputs:* `lem3_1_lcp.csv` (and related experiment data in `lemma31_lcp_powers_experiment.csv`)
* **`lem3_1_lcpPlot.py`**: Reads the LCP CSV and plots the normalized ratio of irreducible LCPs against the theoretical proof ceiling.
* **`reverseBWT.py`**: Verifies Corollary 3.11, which relates the runs of a forward string ($r$) to the runs of its reverse ($\bar{r}$).
    * *Outputs:* `reverse_bwt_experiment.csv`
* **`reverseBWTPlot.py`**: Reads the reverse BWT CSV to verify that $\bar{r} / (r \log_2^2 n) \le O(1)$.

---

## Generated Artifacts

### Data Files (CSVs)
* **`varyingLen.csv`**: Contains `[String_Length, Start_Index, BWT_Runs_r, LZ77_Phrases_z, BWT_String]` for DNA strings of increasing lengths.
* **`sameLen.csv`**: Contains metrics for thousands of DNA strings fixed at length 2000.
* **`spFamily.csv`**: Contains metrics for the generated theoretical worst-case strings for different alphabet sizes.
* **`smBound.csv`**: Contains substring window sizes ($m$), unique substring counts ($|S_m|$), theoretical bounds ($m \cdot z$), and their ratios.
* **`lem3_1_lcp.csv`** / **`lemma31_lcp_powers_experiment.csv`**: Contains counts of irreducible LCPs grouped by bucket ranges $[l, 2l)$.
* **`reverse_bwt_experiment.csv`**: Contains string lengths alongside forward runs ($r$) and reverse runs ($\bar{r}$).

### Visualizations (`res/` directory)
All plots are saved in the `res/` directory.

**From `plots.py`:**
* `r_vs_n.png`: Growth of BWT runs ($r$) as string length ($n$) increases.
* `z_vs_n.png`: Growth of LZ77 phrases ($z$) as string length ($n$) increases.
* `r_vs_z.png`: Direct correlation plot mapping $r$ against $z$.
* `r_over_z_vs_n.png`: How the ratio $r/z$ behaves as $n$ scales.
* `bound_verification.png`: The ultimate check of the upper bound: $r / (z \log_2^2 n)$ across all datasets.
* `bound_verification_alpha4.png`: A cleaner view of the upper bound comparing only DNA against the worst-case (Alphabet=4) string.

**From Specific Lemmas:**
* `sm_bound.png` (via `smBound_plot.py`): Absolute counts of $|S_m|$ vs the theoretical upper bound $m \cdot z$.
* `sm_bound_ratio.png` (via `smBound_plot.py`): Ratio of $|S_m| / (m \cdot z)$ showing it stays below the theoretical maximum of 1.0.
* `lemma31_binned_ratio.png` (via `lem3_1_lcpPlot.py`): Bar chart showing the distribution of irreducible LCPs across powers of 2, compared against the proof's upper bound ceiling of 6.0.
* `lemma31_lcp_counts.png` (via `lem3_1_lcpPlot.py`): Raw counts of irreducible LCPs.
* `lemma31_lcp_ratio.png` (via `lem3_1_lcpPlot.py`): Detailed ratio plotting for Lemma 3.1.
* `corollary311_r_vs_r_rev.png` (via `reverseBWTPlot.py`): Scatter plot mapping forward runs against reverse runs.
* `corollary311_bound.png` (via `reverseBWTPlot.py`): Verification that the normalized ratio between forward and reverse runs remains bounded by $O(1)$ as length increases.

## Execution Order
1. Ensure `y_chromo.fasta` is in the root directory.
2. Run the generation scripts (`varyingLen.py`, `samLen.py`, `spFamily.py`, `smBound.py`, `lem3_1_lcp.py`, `reverseBWT.py`). 
3. Run the plotting scripts (`plots.py`, `smBound_plot.py`, `lem3_1_lcpPlot.py`, `reverseBWTPlot.py`).
4. Check the `res/` folder for your results.