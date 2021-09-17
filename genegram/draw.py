from collections import namedtuple
from pathlib import Path
from typing import List, Dict

from PIL import Image, ImageDraw
from cfpq_data import cfg_from_txt
from tqdm import tqdm

from genegram.cfpq_pyalgo import CNF, BooleanMatrixGraph, all_pairs_reachability_matrix
from genegram.shared import ROOT

__all__ = [
    "seq_fasta_to_pictures",
]

NUCLEOTIDE_TO_COLOR = {"a": 32, "c": 64, "g": 96, "u": 128}

Edge = namedtuple("Edge", ["node_from", "label", "node_to"])


def seq_fasta_to_graphs(seq_fasta: Path) -> Dict[int, List[Edge]]:
    graphs = dict()
    with open(seq_fasta, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                id = line.strip()[1:]
            elif len(id) > 0:
                graphs[int(id)] = [
                    Edge(i, c, i + 1) for i, c in enumerate(line.strip().lower())
                ]
    return graphs


def seq_fasta_to_pictures(seq_fasta: Path, target_dir: Path) -> None:
    # create dir for `target_dir`
    data_dir = Path(target_dir).resolve()
    data_dir.mkdir(parents=True, exist_ok=True)

    # load grammar
    path_to_grammar = str(ROOT / "grammar.txt")
    grammar = CNF.from_cfg(cfg_from_txt(path_to_grammar))

    # translate `seq_fasta` file to graphs
    graphs = seq_fasta_to_graphs(seq_fasta)

    for rna_id, triples in tqdm(graphs.items(), desc="RNA to PNG"):
        bmg = BooleanMatrixGraph.from_triples(triples)

        reachabilities = all_pairs_reachability_matrix(bmg, grammar)

        reachability_pairs = set()
        I, J, V = reachabilities.to_lists()
        for i, _ in enumerate(V):
            reachability_pairs.add((I[i], J[i]))

        # create white&black 8-bit image
        im = Image.new(mode="L", size=(bmg.matrices_size - 1, bmg.matrices_size - 1))
        im_draw = ImageDraw.Draw(im)

        # draw reachabilities
        for i, j in reachability_pairs:
            im_draw.line(xy=[(j - 3, i + 2), (j - 1, i)], fill=255)

        # draw letters
        for i, label, j in triples:
            im_draw.point(xy=(j - 1, i), fill=NUCLEOTIDE_TO_COLOR[label])

            im.save(data_dir / f"{rna_id}.png", "png")
