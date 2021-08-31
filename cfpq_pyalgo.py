from typing import AbstractSet, Iterable, Tuple

from cfpq_data import cnf_from_cfg
from networkx import MultiDiGraph
from pyformlang.cfg import CFG, Variable, Terminal
from pygraphblas import Matrix, BOOL


class CNF:
    def __init__(
        self,
        start_symbol: Variable,
        variables: AbstractSet[Variable],
        terminals: AbstractSet[Terminal],
        unary_productions: Iterable[Tuple[Variable, Terminal]],
        double_productions: Iterable[Tuple[Variable, Variable, Variable]],
    ):
        self.start_symbol = start_symbol
        self.variables = variables
        self.terminals = terminals
        self.unary_productions = unary_productions
        self.double_productions = double_productions

    @classmethod
    def from_cfg(cls, cfg: CFG):
        base_cnf = cnf_from_cfg(cfg)

        unary_productions = list()
        double_productions = list()

        for p in base_cnf.productions:
            if len(p.body) == 0:
                unary_productions.append((p.head, Terminal("$")))
            elif len(p.body) == 1:
                unary_productions.append((p.head, Terminal(p.body[0].value)))
            elif len(p.body) == 2:
                double_productions.append(
                    (p.head, Variable(p.body[0].value), Variable(p.body[1].value))
                )

        cnf = CNF(
            base_cnf.start_symbol,
            base_cnf.variables,
            base_cnf.terminals,
            unary_productions,
            double_productions,
        )

        return cnf

    @classmethod
    def from_text(cls, text, start_symbol: Variable = Variable("S")):
        return CNF.from_cfg(CFG.from_text(text, start_symbol))


class BooleanMatrixGraph:
    def __init__(self, matrices_size: int):
        self.matrices_size = matrices_size
        self.matrices = dict()

    def __getitem__(self, item) -> Matrix:
        if item not in self.matrices:
            self.matrices[item] = Matrix.sparse(
                BOOL, self.matrices_size, self.matrices_size
            )
        return self.matrices[item]

    def __setitem__(self, key, value):
        self.matrices[key] = value

    def __iter__(self):
        return self.matrices.__iter__()

    @property
    def labels(self):
        return list(self.matrices.keys())

    @classmethod
    def from_multidigraph(cls, g: MultiDiGraph):
        bmg = BooleanMatrixGraph(g.number_of_nodes())

        for u, v, edge_labels in g.edges(data=True):
            bmg[Terminal(edge_labels["label"])][u, v] = 1

        return bmg

    @classmethod
    def from_triples(cls, triples):
        number_of_nodes = max({max(u, v) for u, label, v in triples})

        bmg = BooleanMatrixGraph(number_of_nodes + 1)

        for u, label, v in triples:
            bmg[Terminal(label)][u, v] = 1

        return bmg


def all_pairs_reachability_matrix(graph: BooleanMatrixGraph, grammar: CNF):
    m = BooleanMatrixGraph(graph.matrices_size)
    for l, r in grammar.unary_productions:
        m[l] += graph[r]

    changed = True
    while changed:
        changed = False
        for l, r1, r2 in grammar.double_productions:
            old_nnz = m[l].nvals
            m[l] += m[r1].mxm(m[r2], semiring=BOOL.LOR_LAND)
            new_nnz = m[l].nvals

            if old_nnz != new_nnz:
                changed = True

    return m[grammar.start_symbol]
