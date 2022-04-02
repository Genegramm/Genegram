from pyformlang.cfg import CFG

from genegram import BooleanMatrixGraph, WCNF, all_pairs_reachability_matrix


def test_empty_graph():
    bmg = BooleanMatrixGraph()
    wcnf = WCNF.from_text("S -> a")

    I, J, _ = all_pairs_reachability_matrix(bmg, wcnf).to_lists()

    assert set(zip(I, J)) == set()


def test_no_reachable_vertices():
    bmg = BooleanMatrixGraph()
    bmg.add_edge(0, 1, "b")
    bmg.add_edge(1, 2, "c")
    bmg.add_edge(2, 0, "d")

    wcnf = WCNF.from_text("S -> S S | a")

    I, J, _ = all_pairs_reachability_matrix(bmg, wcnf).to_lists()

    assert set(zip(I, J)) == set()


def test_worst_case():
    # The graph is two cycles of coprime lengths with a single common vertex
    # The first cycle is labeled by the open bracket
    # The second cycle is labeled by the close bracket
    bmg = BooleanMatrixGraph()
    bmg.add_edge(0, 1, "a")
    bmg.add_edge(1, 2, "a")
    bmg.add_edge(2, 0, "a")
    bmg.add_edge(0, 3, "b")
    bmg.add_edge(3, 0, "b")

    wcnf = WCNF.from_text("S -> a S b | a b")

    I, J, _ = all_pairs_reachability_matrix(bmg, wcnf).to_lists()

    assert set(zip(I, J)) == {(0, 0), (1, 0), (0, 3), (1, 3), (2, 0), (2, 3)}


def test_full_graph_result():
    # The case when the input graph is sparse, but the result is a full graph.
    # Input graph is a cycle, all edges of which are labeled by the same token
    bmg = BooleanMatrixGraph()
    bmg.add_edge(0, 1, label="a")
    bmg.add_edge(1, 2, label="a")
    bmg.add_edge(2, 3, label="a")
    bmg.add_edge(3, 0, label="a")

    wcnf = WCNF.from_text("S -> S S | a")

    I, J, _ = all_pairs_reachability_matrix(bmg, wcnf).to_lists()

    assert set(zip(I, J)) == {
        (v, to) for v in range(bmg.number_of_nodes) for to in range(bmg.number_of_nodes)
    }
