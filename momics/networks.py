import itertools
import pandas as pd
import networkx as nx
from typing import Dict, List, Tuple


def interaction_to_graph(
    df: pd.DataFrame, pos_cutoff: float = 0.8, neg_cutoff: float = -0.6
) -> Tuple[List[str], List[Tuple[str, str]], List[Tuple[str, str]]]:
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


def pairwise_jaccard_lower_triangle(
    network_results: dict, edge_type: str = "all"
) -> pd.DataFrame:
    """
    Calculate pairwise Jaccard similarity for the lower triangle of all group comparisons.
    Returns a DataFrame with columns: group1, group2, jaccard_similarity.

    If `edge_type` is 'all', it calculates Jaccard similarity for all edges in the graphs.

    Args:
        network_results (dict): Dictionary containing network results for each group.
            Keys are 'graph', 'nodes', and lists of specific edges from the graph.
        edge_type (str): dict key for list of edges to consider (or 'all').

    Returns:
        pd.DataFrame: DataFrame containing pairwise Jaccard similarity.
    """
    # Extract all group names
    groups = list(network_results.keys())
    results = []

    # define empty DataFrame with groups as index and columns
    pivoted = pd.DataFrame(index=groups, columns=groups)
    # Iterate over all unique pairs (lower triangle, i < j)
    for g1, g2 in itertools.combinations(groups, 2):
        if edge_type == "all":
            edges1 = set(network_results[g1]["graph"].edges())
            edges2 = set(network_results[g2]["graph"].edges())
        else:
            edges1 = set(network_results[g1][edge_type])
            edges2 = set(network_results[g2][edge_type])
        intersection = edges1 & edges2
        union = edges1 | edges2
        jaccard = len(intersection) / len(union) if len(union) > 0 else float("nan")
        results.append({"group1": g1, "group2": g2, "jaccard_similarity": jaccard})
        pivoted.loc[g2, g1] = jaccard

    return pivoted


def build_interaction_graphs(
    correlation_data: dict,
    pos_cutoff: float = 0.5,
    neg_cutoff: float = -0.5,
    p_val_cutoff: float = 0.05,
) -> Dict:
    """
    Build interaction graphs from correlation data.

    Args:
        correlation_data (dict): A dictionary containing correlation data for different factors.
        pos_cutoff (float): The positive correlation cutoff.
        neg_cutoff (float): The negative correlation cutoff.
        p_val_cutoff (float): The p-value cutoff.

    Returns:
        Dict: A dictionary containing network results for each factor.
    """
    network_results = {}
    for factor, dict_df in correlation_data.items():
        print(f"Factor: {factor}")
        nodes, edges_pos, edges_neg = interaction_to_graph_with_pvals(
            dict_df["correlation"],
            dict_df["p_vals_fdr"],
            pos_cutoff=pos_cutoff,
            neg_cutoff=neg_cutoff,
            p_val_cutoff=p_val_cutoff,
        )
        G = nx.Graph(mode=factor)

        G.add_nodes_from(nodes)
        G.add_edges_from(edges_pos, color="green")
        G.add_edges_from(edges_neg, color="red")

        network_results[factor] = {
            "graph": G,
            "nodes": nodes,
            "edges_pos": edges_pos,
            "edges_neg": edges_neg,
        }

        degree_centrality = nx.degree_centrality(G)

        network_results[factor]["degree_centrality"] = sorted(
            degree_centrality.items(), key=lambda x: x[1], reverse=True
        )[:10]

        betweenness = nx.betweenness_centrality(G)

        network_results[factor]["top_betweenness"] = sorted(
            betweenness.items(), key=lambda x: x[1], reverse=True
        )[:10]
        network_results[factor]["bottom_betweenness"] = sorted(
            betweenness.items(), key=lambda x: x[1]
        )[:10]
        network_results[factor]["total_nodes"] = G.number_of_nodes()
        network_results[factor]["total_edges"] = G.number_of_edges()

    return network_results
