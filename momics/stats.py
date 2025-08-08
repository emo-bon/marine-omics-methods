import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict
from scipy.stats import spearmanr, pearsonr


def spearman_from_taxonomy(split_taxonomy: Dict) -> Dict:
    """
    Compute Spearman correlation and p-values for the full taxonomy split by a factor.
    Refer `momics.taxonomy.split_taxonomic_data` for more information.

    Args:
        split_taxonomy (dict): A dictionary containing dataframes for each factor.

    Returns:
        dict: A dictionary containing Spearman correlation and p-values for each factor.
    """
    spearman_taxa = {}
    # Compute Spearman correlation
    for factor, df in split_taxonomy.items():
        corr, p_spearman = spearmanr(df.T)
        assert (
            corr.shape == p_spearman.shape
        ), "Spearman correlation and p-values must have the same shape."

        spearman_taxa[factor] = {
            "correlation": pd.DataFrame(corr, index=df.index, columns=df.index),
            "p_vals": pd.DataFrame(p_spearman, index=df.index, columns=df.index),
        }
        assert (
            spearman_taxa[factor]["correlation"].shape
            == spearman_taxa[factor]["p_vals"].shape
        ), "Spearman correlation and p-values must have the same shape."

    return spearman_taxa


################################
## Plotting correlations etc. ##
################################
def plot_association_histogram(assoc_data: Dict, bins: int = 50) -> None:
    """
        Plot a histogram of the correlation values for each factor.

    Args:
        correlation_data (dict): A dictionary containing correlation data for each factor.
        bins (int): The number of bins to use for the histogram.

    Returns:
        None
    """
    # histogram of the correlation values for setting graph cutoffs
    plt.figure(figsize=(10, 5))
    for factor, df in assoc_data.items():
        try:  # for bray-curtis
            lower_triangle = df.values[np.tril_indices_from(df.values, k=-1)]
        except Exception as e:
            lower_triangle = df["correlation"].values[
                np.tril_indices_from(df["correlation"].values, k=-1)
            ]

        plt.hist(lower_triangle, bins=bins, alpha=0.5, label=factor, log=True)
        print("values in the triangle:", len(lower_triangle))

    plt.title("Histogram of Correlation Values")
    plt.xlabel("Correlation")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()


def plot_fdr(correlations: Dict, pval_cutoff: float) -> None:
    """
    Plot the FDR-corrected p-values against the raw p-values for each factor.

    Args:
        correlations (dict): A dictionary containing correlation data for each factor.
        pval_cutoff (float): The p-value cutoff for significance.

    Returns:
        None
    """
    plt.figure(figsize=(10, 5))
    for factor, df in correlations.items():
        raw_pvals = df["p_vals"]
        pvals = df["p_vals_fdr"]
        significant_raw = raw_pvals < pval_cutoff
        significant_fdr = pvals < pval_cutoff

        plt.scatter(raw_pvals, pvals, alpha=0.7, label=factor)

        # Optionally, print how many associations were removed by FDR
        removed = np.sum(
            significant_raw.values.flatten() & ~significant_fdr.values.flatten()
        )
        significant = np.sum(significant_fdr.values.flatten())
        print(f"### FACTOR: {factor} ###")
        print(
            f"Significant associations before {np.sum(significant_raw.values.flatten())} and after {significant}"
        )
        print(f"Associations significant before FDR but not after: {removed}")

    plt.axvline(pval_cutoff, color="gray", linestyle="--", label=f"Raw p={pval_cutoff}")
    plt.axhline(
        pval_cutoff, color="black", linestyle="--", label=f"FDR p={pval_cutoff}"
    )
    plt.xlabel("Raw p-value")
    plt.ylabel("FDR-corrected p-value")
    plt.title(f"Comparison of Raw and FDR-corrected p-values for {factor}")
    plt.legend()
    plt.show()
