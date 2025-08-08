import pytest
import pandas as pd
import networkx as nx
import os

from momics.metadata import (
    get_metadata_udal,
    enhance_metadata,
    clean_metadata,
)
from momics.networks import (
    interaction_to_graph,
    interaction_to_graph_with_pvals,
    pairwise_jaccard_lower_triangle,
)
from momics.constants import COL_NAMES_HASH_EMO_BON_VRE as COL_NAMES_HASH
from momics.constants import TAXONOMY_RANKS


def test_interaction_to_graph():
    # Create a 3x3 correlation matrix
    corr = pd.DataFrame(
        [
            [1.0, 0.9, -0.7],
            [0.9, 1.0, 0.2],
            [-0.7, 0.2, 1.0],
        ],
        columns=["A", "B", "C"],
        index=["A", "B", "C"],
    )

    nodes, edges_pos, edges_neg = interaction_to_graph(
        corr, pos_cutoff=0.8, neg_cutoff=-0.6
    )

    # Should find one positive edge (A,B) and one negative edge (A,C)
    assert set(nodes) == {"A", "B", "C"}
    assert ("A", "B") in edges_pos or ("B", "A") in edges_pos
    assert ("A", "C") in edges_neg or ("C", "A") in edges_neg
    assert len(edges_pos) == 1
    assert len(edges_neg) == 1

    # If cutoffs are extreme, no edges should be found
    nodes2, edges_pos2, edges_neg2 = interaction_to_graph(
        corr, pos_cutoff=1.1, neg_cutoff=-1.1
    )
    assert len(edges_pos2) == 0
    assert len(edges_neg2) == 0


def test_interaction_to_graph_with_pvals():
    # Create a 3x3 correlation matrix and matching p-value matrix
    corr = pd.DataFrame(
        [
            [1.0, 0.9, -0.7],
            [0.9, 1.0, 0.2],
            [-0.7, 0.2, 1.0],
        ],
        columns=["A", "B", "C"],
        index=["A", "B", "C"],
    )
    pvals = pd.DataFrame(
        [
            [0.01, 0.01, 0.01],
            [0.01, 0.01, 0.01],
            [0.01, 0.01, 0.01],
        ],
        columns=["A", "B", "C"],
        index=["A", "B", "C"],
    )

    nodes, edges_pos, edges_neg = interaction_to_graph_with_pvals(
        corr, pvals, pos_cutoff=0.8, neg_cutoff=-0.6, p_val_cutoff=0.05
    )

    # Should find one positive edge (A,B) and one negative edge (A,C)
    assert set(nodes) == {"A", "B", "C"}
    assert ("A", "B") in edges_pos or ("B", "A") in edges_pos
    assert ("A", "C") in edges_neg or ("C", "A") in edges_neg
    assert len(edges_pos) == 1
    assert len(edges_neg) == 1

    # If p-value is too high, no edges should be found
    pvals_high = pvals + 1
    _, edges_pos2, edges_neg2 = interaction_to_graph_with_pvals(
        corr, pvals_high, pos_cutoff=0.8, neg_cutoff=-0.6, p_val_cutoff=0.05
    )
    assert len(edges_pos2) == 0
    assert len(edges_neg2) == 0


def test_pairwise_jaccard_lower_triangle():
    # Create mock graphs for three groups
    g1 = nx.Graph()
    g1.add_edges_from([("A", "B"), ("B", "C")])
    g2 = nx.Graph()
    g2.add_edges_from([("A", "B"), ("C", "D")])
    g3 = nx.Graph()
    g3.add_edges_from([("X", "Y")])

    network_results = {
        "group1": {"graph": g1},
        "group2": {"graph": g2},
        "group3": {"graph": g3},
    }

    df = pairwise_jaccard_lower_triangle(network_results, edge_type="all")

    # Check DataFrame shape and index/columns
    assert isinstance(df, pd.DataFrame)
    assert set(df.index) == {"group1", "group2", "group3"}
    assert set(df.columns) == {"group1", "group2", "group3"}

    # Jaccard between group1 and group2: 1 shared edge, 3 total edges
    assert df.loc["group2", "group1"] == pytest.approx(1 / 3)
    # Jaccard between group1 and group3: 0 shared, 3 total
    assert df.loc["group3", "group1"] == 0.0
    # Jaccard between group2 and group3: 0 shared, 3 total
    assert df.loc["group3", "group2"] == 0.0
    # Upper triangle and diagonal should be None/NaN
    assert pd.isna(df.loc["group1", "group2"])
    assert pd.isna(df.loc["group1", "group1"])
