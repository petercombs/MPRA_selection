from argparse import ArgumentParser
from dendropy import Tree
from Bio import SeqIO, AlignIO
import pandas as pd
from os import path


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--target-species-fasta",
        "-t",
        default=False,
        help="""Fasta file containing only the species we're interested in (for
        example, only primates)""",
    )
    parser.add_argument("--enhancer-name", "-e", default=None)
    parser.add_argument("seqs")
    parser.add_argument("input_tree")
    parser.add_argument("mpra_data")
    args = parser.parse_args()
    if args.enhancer_name is None:
        # By default, we're going to assume that the ancestral reconstruction
        # program has its own directory within the enhancer.
        args.enhancer_name = path.basename(path.dirname(path.dirname(args.input_tree)))
    return args


if __name__ == "__main__":
    args = parse_args()
    tree = Tree.get_from_path(args.input_tree, "newick")
    alignment = AlignIO.read(args.seqs, "fasta")
    alignment_posns = {rec.id: i for i, rec in enumerate(alignment)}

    mpra_data = pd.read_csv(
        args.mpra_data,
        sep="\t",
        header=None,
        names=["Enhancer", "pos", "newbase", "effect", "pval"],
        index_col=["Enhancer", "pos", "newbase"],
    )

    print(tree.as_string("newick"))
    if args.target_species_fasta:
        target_species = {
            rec.id.replace("_", " ")
            for rec in SeqIO.parse(args.target_species_fasta, "fasta")
        }
        tree.retain_taxa_with_labels(target_species)

    print(tree.as_string("newick"))
    for node in tree:
        if node.parent_node is None:
            continue
        parent_name = node.parent_node.label.replace(" ", "_")
        node_name = (node.label or node.taxon.label).replace(" ", "_")

        align_len = len(alignment[0])
        homo_pos = 0
        homo_seq = alignment[alignment_posns["Homo_sapiens"]]
        parent_seq = alignment[alignment_posns[parent_name]]
        child_seq = alignment[alignment_posns[node_name]]
        du = 0
        dd = 0
        dn = 0
        for i in range(align_len):
            if homo_seq[i] == "-":
                continue
            homo_pos += 1
            if parent_seq[i] == "-" or child_seq[i] == "-":
                continue
            elif parent_seq[i] == homo_seq[i] and parent_seq[i] != child_seq[i]:
                mpra_row = mpra_data.loc[args.enhancer_name, homo_pos, child_seq[i]]
                if mpra_row.pval > .05:
                    dn += 1
                elif mpra_row.effect < 0:
                    dd += 1
                else:
                    du += 1
            elif child_seq[i] == homo_seq[i] and parent_seq[i] != child_seq[i]:
                mpra_row = mpra_data.loc[args.enhancer_name, homo_pos, parent_seq[i]]
                if mpra_row.pval > .05:
                    dn += 1
                elif mpra_row.effect > 0:
                    dd += 1
                else:
                    du += 1

        print(parent_name, node_name, du, dn, '<-- UP')
        print(parent_name, node_name, dd, dn, '<-- DOWN')
