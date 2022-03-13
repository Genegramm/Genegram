from pyformlang.cfg import CFG, Production, Variable, Terminal

from genegram import WCNF


def test_empty():
    cfg = CFG.from_text("S -> epsilon")

    wcnf = WCNF(cfg)

    assert wcnf.start_variable == Variable("S")
    assert wcnf.variables == [Variable("S")]
    assert wcnf.terminals == []
    assert wcnf.productions == [Production(Variable("S"), [])]
    assert wcnf.epsilon_productions == [Production(Variable("S"), [])]
    assert wcnf.unary_productions == []
    assert wcnf.binary_productions == []


def test_a():
    cfg = CFG.from_text("S -> a")

    wcnf = WCNF(cfg)

    assert wcnf.start_variable == Variable("S")
    assert wcnf.variables == [Variable("S")]
    assert wcnf.terminals == [Terminal("a")]
    assert wcnf.productions == [Production(Variable("S"), [Terminal("a")])]
    assert wcnf.epsilon_productions == []
    assert wcnf.unary_productions == [Production(Variable("S"), [Terminal("a")])]
    assert wcnf.binary_productions == []
