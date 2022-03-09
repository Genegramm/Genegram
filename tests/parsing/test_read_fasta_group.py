from pathlib import Path

from genegram import RNA, read_fasta_group

root = Path(__file__).parent.resolve()
fasta = root.parent / "data" / "seq.fasta"


def test_limit_0():
    data = list(read_fasta_group(fasta, 0))

    assert data == [
        [RNA(description="34551", sequence="GGCCUCCAAGCUGUGCCUUGGGUGGCC")],
        [RNA(description="34552", sequence="CCUCCCUUACAAGGAGG")],
        [RNA(description="34553", sequence="GGAGUGGCCGAAAGGCAUCUCC")],
        [RNA(description="34735", sequence="GGCUCUCAGUGAGCC")],
        [RNA(description="34736", sequence="GGUGCUCAGUAGGAGACGAACCGCACC")],
        [RNA(description="34737", sequence="GGGCAAGCCC")],
        [RNA(description="34738", sequence="GGCAGUGUGAGUACCUUCACACGUC")],
        [RNA(description="34745", sequence="GGUGUGAACACC")],
        [RNA(description="34767", sequence="GGUGCAUGGCACC")],
        [
            RNA(
                description="34840",
                sequence="GGCGGUACUAGUUGAGAAACUAGCUCUGUAUCUGGCGGACCCGUGGUGGAACUGUGAAGUUCGGAACACCCGGCCGCAACCCUGGGAGAGGUCCCAGGGUU",
            )
        ],
    ]


def test_limit_50():
    data = list(read_fasta_group(fasta, 50))

    assert data == [
        [
            RNA(description="34551", sequence="GGCCUCCAAGCUGUGCCUUGGGUGGCC"),
            RNA(description="34552", sequence="CCUCCCUUACAAGGAGG"),
        ],
        [
            RNA(description="34553", sequence="GGAGUGGCCGAAAGGCAUCUCC"),
            RNA(description="34735", sequence="GGCUCUCAGUGAGCC"),
        ],
        [
            RNA(description="34736", sequence="GGUGCUCAGUAGGAGACGAACCGCACC"),
            RNA(description="34737", sequence="GGGCAAGCCC"),
        ],
        [
            RNA(description="34738", sequence="GGCAGUGUGAGUACCUUCACACGUC"),
            RNA(description="34745", sequence="GGUGUGAACACC"),
            RNA(description="34767", sequence="GGUGCAUGGCACC"),
        ],
        [
            RNA(
                description="34840",
                sequence="GGCGGUACUAGUUGAGAAACUAGCUCUGUAUCUGGCGGACCCGUGGUGGAACUGUGAAGUUCGGAACACCCGGCCGCAACCCUGGGAGAGGUCCCAGGGUU",
            )
        ],
    ]


def test_limit_500():
    data = list(read_fasta_group(fasta, 500))

    assert data == [
        [
            RNA(description="34551", sequence="GGCCUCCAAGCUGUGCCUUGGGUGGCC"),
            RNA(description="34552", sequence="CCUCCCUUACAAGGAGG"),
            RNA(description="34553", sequence="GGAGUGGCCGAAAGGCAUCUCC"),
            RNA(description="34735", sequence="GGCUCUCAGUGAGCC"),
            RNA(description="34736", sequence="GGUGCUCAGUAGGAGACGAACCGCACC"),
            RNA(description="34737", sequence="GGGCAAGCCC"),
            RNA(description="34738", sequence="GGCAGUGUGAGUACCUUCACACGUC"),
            RNA(description="34745", sequence="GGUGUGAACACC"),
            RNA(description="34767", sequence="GGUGCAUGGCACC"),
            RNA(
                description="34840",
                sequence="GGCGGUACUAGUUGAGAAACUAGCUCUGUAUCUGGCGGACCCGUGGUGGAACUGUGAAGUUCGGAACACCCGGCCGCAACCCUGGGAGAGGUCCCAGGGUU",
            ),
        ]
    ]
