from collections import defaultdict
from typing import Dict, List

import pygraphblas as gb
from pyformlang.cfg import CFG, Production, Variable, Terminal, Epsilon

__all__ = [
    "all_pairs_reachability_matrix",
    "WCNF",
    "BooleanMatrixGraph",
]


class BooleanMatrixGraph:
    def __init__(self, matrices_size):
        self._matrices: Dict[str, gb.Matrix] = dict()
        self._matrices_size: int = matrices_size

    def __getitem__(self, label) -> gb.Matrix:
        if label not in self._matrices:
            self._matrices[label] = gb.Matrix.sparse(
                typ=gb.BOOL,
                nrows=self._matrices_size,
                ncols=self._matrices_size,
            )
        return self._matrices[label]

    def __setitem__(self, label, matrix: gb.Matrix) -> None:
        self._matrices[label] = matrix


class WCNF:
    def __init__(self, cfg: CFG):
        self._cfg: CFG = cfg
        self.start_symbol: Variable = cfg.start_symbol

        if not _is_in_wcnf(cfg):
            cnf = cfg.to_normal_form()
        else:
            cnf = cfg

        self.epsilon_productions: List[Production] = []
        self.unary_productions: List[Production] = []
        self.binary_productions: List[Production] = []

        for production in self._cfg.productions:
            if production.body in ([], Epsilon):
                if production not in self.epsilon_productions:
                    self.epsilon_productions.append(production)

        for production in cnf.productions:
            if len(production.body) == 1:
                if production not in self.unary_productions:
                    self.unary_productions.append(production)
            elif len(production.body) == 2:
                if production not in self.binary_productions:
                    self.binary_productions.append(production)

        self.productions = (
            self.epsilon_productions + self.unary_productions + self.binary_productions
        )

        self.variables: List[Variable] = []
        self.terminals: List[Terminal] = []

        for production in self.productions:
            if production.head not in self.variables:
                self.variables.append(production.head)

            for term in production.body:
                if isinstance(term, Terminal):
                    if term not in self.terminals:
                        self.terminals.append(term)
                elif isinstance(term, Variable):
                    if term not in self.variables:
                        self.variables.append(term)

        self.variables.sort(key=str)
        self.terminals.sort(key=str)
        self.epsilon_productions.sort(key=str)
        self.unary_productions.sort(key=str)
        self.unary_productions.sort(key=str)
        self.binary_productions.sort(key=str)
        self.productions.sort(key=str)

    def contains(self, word: str):
        return self._cfg.contains(word)


def _is_in_wcnf(cfg: CFG) -> bool:
    for production in cfg.productions:
        if len(production.body) > 2:
            return False
        elif len(production.body) == 2:
            if not (
                isinstance(production.body[0], Variable)
                and isinstance(production.body[1], Variable)
            ):
                return False
        elif len(production.body) == 1:
            if not isinstance(production.body[0], Terminal):
                return False
    return True


def all_pairs_reachability_matrix(graph: BooleanMatrixGraph, grammar: WCNF):
    m = BooleanMatrixGraph(graph._matrices_size)

    for production in grammar.unary_productions:
        # production :: l -> r
        l = production.head
        r = production.body[0]

        m[l] += graph[r]

    changed = True
    nvals = defaultdict(int)
    while changed:
        changed = False
        for production in grammar.binary_productions:
            # production :: l -> r1 r2
            l = production.head
            r1 = production.body[0]
            r2 = production.body[1]

            m[l] += m[r1].mxm(m[r2], semiring=gb.BOOL.ANY_PAIR)

            nnz = m[l].nvals

            changed |= nvals[l] != nnz

            nvals[l] = nnz

    return m[grammar.start_symbol]
