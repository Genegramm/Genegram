from argparse import ArgumentParser
from pathlib import Path

from PIL import Image, ImageDraw
from cfpq_data import cfg_from_txt
from tqdm import tqdm

from cfpq_pyalgo import CNF, BooleanMatrixGraph, all_pairs_reachability_matrix

NUCLEOTIDE_TO_COLOR = {"a": 32, "c": 64, "g": 96, "u": 128}
ROOT_DIR = Path(__file__).parent.resolve()


def seq_fasta_to_graphs(seq_fasta):
    graphs = dict()
    with open(seq_fasta, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                id = line.strip()[1:]
            elif len(id) > 0:
                graphs[int(id)] = [
                    (i, c, i + 1) for i, c in enumerate(line.strip().lower())
                ]
    return graphs


def seq_fasta_to_pictures(seq_fasta, target_dir):
    # create dir for `target_dir`
    data_dir = Path(target_dir).resolve()
    data_dir.mkdir(parents=True, exist_ok=True)

    # load grammar
    path_to_grammar = str(ROOT_DIR / "grammar.txt")
    grammar = CNF.from_cfg(cfg_from_txt(path_to_grammar))

    # translate `seq_fasta` file to graphs
    graphs = seq_fasta_to_graphs(seq_fasta)

    for rna_id, triples in tqdm(graphs.items()):
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


if __name__ == "__main__":
    parser = ArgumentParser(description="seq.fasta to pictures")
    parser.add_argument("--path", required=True, type=str, help="path to seq.fasta")
    parser.add_argument(
        "--target_dir",
        required=True,
        type=str,
        help="path to the folder in which to save the pictures",
    )

    args = parser.parse_args()

    seq_fasta_to_pictures(args.path, args.target_dir)
