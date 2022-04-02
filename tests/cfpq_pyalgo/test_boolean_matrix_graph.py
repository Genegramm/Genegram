from pygraphblas import Matrix, BOOL

from genegram import BooleanMatrixGraph


def test_empty():
    bmg = BooleanMatrixGraph()

    assert bmg.number_of_nodes == 0
    assert bmg._matrices == dict()


def test_one_edge():
    bmg = BooleanMatrixGraph()
    bmg.add_edge(0, 1, "label")

    assert bmg.number_of_nodes == 2
    assert bmg._matrices == {
        "label": Matrix.from_lists(
            I=[0],
            J=[1],
            V=[True],
            nrows=2,
            ncols=2,
            typ=BOOL,
        )
    }


def test_two_edges():
    bmg = BooleanMatrixGraph()
    bmg.add_edge(0, 1, "A")
    bmg.add_edge(1, 2, "B")

    assert bmg._number_of_nodes == 3
    assert bmg._matrices == {
        "A": Matrix.from_lists(
            I=[0],
            J=[1],
            V=[True],
            nrows=3,
            ncols=3,
            typ=BOOL,
        ),
        "B": Matrix.from_lists(
            I=[1],
            J=[2],
            V=[True],
            nrows=3,
            ncols=3,
            typ=BOOL,
        ),
    }
