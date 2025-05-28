import pandas as pd


def interaction_to_graph(
    df: pd.DataFrame, pos_cutoff: float = 0.8, neg_cutoff: float = -0.6
) -> tuple:
    """
    Create a network from the correlation matrix.
    Args:
        df (pd.DataFrame): The input DataFrame containing correlation values.
        pos_cutoff (float): Positive correlation cutoff.
        neg_cutoff (float): Negative correlation cutoff.
    Returns:
        nodes (list): List of node indices.
        edges_pos (list): List of positive edges.
        edges_neg (list): List of negative edges.
    """
    nodes, edges_pos, edges_neg = [], [], []
    count_pos, count_neg = 0, 0
    cols = df.columns.tolist()
    for i in range(df.shape[0]):
        nodes.append(cols[i])
        for j in range(i + 1, df.shape[1]):
            if df.iloc[i, j] > pos_cutoff:
                edges_pos.append((cols[i], cols[j]))
                count_pos += 1
            # print(f"Sample {i} and Sample {j} have a high correlation of {df.iloc[i, j]}")
            elif df.iloc[i, j] < neg_cutoff:
                edges_neg.append((cols[i], cols[j]))
                count_neg += 1
                # print(f"Sample {i} and Sample {j} have a high negative correlation of {df.iloc[i, j]}")
    print(f"Number of positive edges: {count_pos}")
    print(f"Number of negative edges: {count_neg}")
    return nodes, edges_pos, edges_neg


def interaction_to_graph_with_pvals(
    df: pd.DataFrame,
    pvals_df: pd.DataFrame,
    pos_cutoff: float = 0.8,
    neg_cutoff: float = -0.6,
    p_val_cutoff: float = 0.05,
) -> tuple:
    """
    Create a network from the correlation matrix and p-values.
    Args:
        df (pd.DataFrame): The input DataFrame containing correlation values.
        pvals_df (pd.DataFrame): The DataFrame containing p-values.
        pos_cutoff (float): Positive correlation cutoff.
        neg_cutoff (float): Negative correlation cutoff.
    Returns:
        nodes (list): List of node indices.
        edges_pos (list): List of positive edges with p-values.
        edges_neg (list): List of negative edges with p-values.
    """
    nodes, edges_pos, edges_neg = [], [], []
    count_pos, count_neg = 0, 0
    cols = df.columns.tolist()
    for i in range(df.shape[0]):
        nodes.append(cols[i])
        for j in range(i + 1, df.shape[1]):
            if df.iloc[i, j] > pos_cutoff and pvals_df.iloc[i, j] < p_val_cutoff:
                edges_pos.append((cols[i], cols[j]))
                count_pos += 1
            elif df.iloc[i, j] < neg_cutoff and pvals_df.iloc[i, j] < p_val_cutoff:
                edges_neg.append((cols[i], cols[j]))
                count_neg += 1
    print(f"Number of positive edges: {count_pos}")
    print(f"Number of negative edges: {count_neg}")
    return nodes, edges_pos, edges_neg
