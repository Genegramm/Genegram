"""The All-Pairs CFL-reachability module"""
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
    """A Labeled Graph decomposed into Boolean Matrices"""

    def __init__(self, number_of_nodes: int = 0):
        self._matrices: Dict[str, gb.Matrix] = dict()
        self._number_of_nodes: int = number_of_nodes

    def __getitem__(self, label: str) -> gb.Matrix:
        if label not in self._matrices:
            self._matrices[label] = gb.Matrix.sparse(
                typ=gb.BOOL,
                nrows=self._number_of_nodes,
                ncols=self._number_of_nodes,
            )
        return self._matrices[label]

    def __setitem__(self, label: str, matrix: gb.Matrix) -> None:
        self._matrices[label] = matrix

    @property
    def number_of_nodes(self) -> int:
        """The number of nodes in the graph
        Returns
        -------
        number_of_nodes: int
            Number of nodes in the graph
        """
        return self._number_of_nodes

    def add_edge(self, u: int, v: int, label: str) -> None:
        """Add an edge between `u` and `v` with label `label`.
        The nodes `u` and `v` will be automatically added if they are
        not already in the graph.
        Parameters
        ----------
        u: int
            The tail of the edge
        v: int
            The head of the edge
        label: str
            The label of the edge
        """
        if max(u, v) > self._number_of_nodes - 1:
            self._number_of_nodes = max(u, v) + 1
            for key in self._matrices:
                self._matrices[key].resize(self._number_of_nodes, self._number_of_nodes)

        if label not in self._matrices:
            self._matrices[label] = gb.Matrix.sparse(
                typ=gb.BOOL, nrows=self._number_of_nodes, ncols=self._number_of_nodes
            )

        self._matrices[label][u, v] = True


class WCNF:
    """A Weak Chomsky Normal Form of Context-Free Grammar
    in which products take the following form:
    - A -> B C
    - A -> a
    - A -> epsilon
    where `A`, `B` and `C` are variables; `a` is an arbitrary terminal
    Also known as Weak Chomsky Normal Form

    Parameters
    ----------
    cfg: CFG
        Context-Free Grammar
    """

    def __init__(self, cfg: CFG):
        self._cfg: CFG = cfg
        self.start_variable: Variable = cfg.start_symbol

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

    def contains(self, word: str) -> bool:
        """Gives the membership of a word to the grammar

        Parameters
        ----------
        word : str
            The word to check

        Returns
        ----------
        contains : bool
            Whether word if in the grammar's language or not
        """
        return self._cfg.contains(word)

    @classmethod
    def from_text(cls, text, start_symbol=Variable("S")):
        """
        Read a Weak Chomsky Normal Form Context-Free Grammar from a text.
        The text contains one rule per line.
        The structure of a production is:
        head -> body1 | body2 | ... | bodyn
        where | separates the bodies.
        A variable (or non terminal) begins by a capital letter.
        A terminal begins by a non-capital character
        Terminals and Variables are separated by spaces.
        An epsilon symbol can be represented by epsilon, $, ε, ϵ or Є.
        If you want to have a variable name starting with a non-capital \
        letter or a terminal starting with a capital letter, you can \
        explicitly give the type of your symbol with "VAR:yourVariableName" \
        or "TER:yourTerminalName" (with the quotation marks). For example:
        S -> "TER:John" "VAR:d" a b

        Parameters
        ----------
        text : str
            The text of transform
        start_symbol : str, optional
            The start symbol, S by default

        Returns
        -------
        wcnf : WCNF
            A Weak Chomsky Normal Form Context-Free Grammar
        """
        return cls(CFG.from_text(text, start_symbol))


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


def all_pairs_reachability_matrix(
    graph: BooleanMatrixGraph, grammar: WCNF
) -> gb.Matrix:
    """Determines the pairs of vertices (`u`, `v`)
    where there exists a path from `u` to `v`
    in `graph` and its word is in the language of `grammar`

    Parameters
    ----------
    graph: BooleanMatrixGraph
    grammar: WCNF

    Returns
    -------
    pairs: pygraphblas.Matrix
        Reachability matrix
    """
    bmg = BooleanMatrixGraph(graph.number_of_nodes)

    for production in grammar.unary_productions:
        # production :: l -> r
        l = production.head.value
        r = production.body[0].value

        bmg[l] += graph[r]

    changed = True
    nvals = defaultdict(int)
    while changed:
        changed = False
        for production in grammar.binary_productions:
            # production :: l -> r1 r2
            l = production.head.value
            r1 = production.body[0].value
            r2 = production.body[1].value

            bmg[l] += bmg[r1].mxm(bmg[r2], semiring=gb.BOOL.ANY_PAIR)

            nnz = bmg[l].nvals

            changed |= nvals[l] != nnz

            nvals[l] = nnz

    return bmg[grammar.start_variable.value]
