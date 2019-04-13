from argparse import ArgumentParser
from dendropy import Tree, Taxon
from Bio import SeqIO, AlignIO
import pandas as pd
from os import path
from scipy.stats import fisher_exact
from sys import stderr

tree_labels = ["SKIP"]

for i in range(1, 10):
    for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        tree_labels.append(c * i)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--target-species-fasta",
        "-t",
        default=False,
        help="""Fasta file containing only the species we're interested in (for
        example, only primates)""",
    )
    parser.add_argument("--output-tree", "-O", default=False)
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
    AlignIO.write(
        alignment, path.join(path.dirname(args.seqs), "aln.clustal"), "clustal"
    )
    alignment_posns = {rec.id: i for i, rec in enumerate(alignment)}

    if "GRCh38_ALL" in args.mpra_data:
        mpra_data = pd.read_csv(args.mpra_data, sep="\t")
        old_mpra_data = mpra_data

        # We are interested in the position per-enhancer
        mpra_data["pos"] = -1
        for enh in mpra_data.Element.unique():
            ix = mpra_data.Element == enh
            mpra_data.loc[ix, "pos"] = (
                mpra_data.loc[ix].Position - mpra_data.loc[ix].Position.min() + 1
            )

        mpra_data = pd.DataFrame(
            index=pd.MultiIndex.from_arrays(
                [mpra_data.Element, mpra_data.pos, mpra_data.Alt]
            ),
            data={"Value": list(mpra_data.Value), "pval": list(mpra_data["P-Value"])},
            # We need to cast to a list to get it to work with the reindexing
        )
    else:
        # This is likely the Patwardhan data
        mpra_data = pd.read_csv(
            args.mpra_data,
            sep="\t",
            header=None,
            names=["Element", "pos", "Alt", "Value", "pval"],
            index_col=["Element", "pos", "Alt"],
        )

    print(tree.as_string("newick"))
    if args.target_species_fasta:
        target_species = {
            rec.id.replace("_", " ")
            for rec in SeqIO.parse(args.target_species_fasta, "fasta")
        }
        outgroup_tree = tree.extract_tree_without_taxa_labels(target_species)
        outgroup_root = tree.find_node(
            lambda n: n.label == outgroup_tree.nodes()[0].label
        )
        outgroup_root.clear_child_nodes()
        outgroup_root.taxon = Taxon(outgroup_root.label)
        tree.taxon_namespace.append(outgroup_root.taxon)

    i = 0
    for edge in tree.edges():
        if edge.head_node is outgroup_root:
            edge.label = "N/A"
            continue
        edge.label = tree_labels[i]
        i += 1

    for node in tree:
        node.annotations["comparison"] = node.incident_edges()[-1].label
        if args.target_species_fasta and node is outgroup_root:
            node.annotations["species"] = 'Outgroup'
        else:
            node.annotations["species"] = node.taxon.label if node.taxon else node.label

    if args.output_tree:
        tree.write_to_path(args.output_tree, "nexus", suppress_annotations=False)

    overall_du = 0
    overall_dd = 0
    overall_dn = 0

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
            elif homo_seq[i] == "N":
                homo_pos += 1
                continue
            homo_pos += 1
            try:
                if parent_seq[i] == "-" or child_seq[i] == "-":
                    continue
                elif parent_seq[i] == homo_seq[i] and parent_seq[i] != child_seq[i]:
                    if len(mpra_data.loc[args.enhancer_name, homo_pos]) < 3:
                        continue
                    mpra_row = mpra_data.loc[args.enhancer_name, homo_pos, child_seq[i]]
                    if mpra_row.pval > .05:
                        dn += 1
                    elif mpra_row.Value < 0:
                        dd += 1
                    else:
                        du += 1
                elif child_seq[i] == homo_seq[i] and parent_seq[i] != child_seq[i]:
                    if len(mpra_data.loc[args.enhancer_name, homo_pos]) < 3:
                        raise ValueError(f"Missing one or more bases at {homo_pos}")
                    mpra_row = mpra_data.loc[
                        args.enhancer_name, homo_pos, parent_seq[i]
                    ]
                    if mpra_row.pval > .05:
                        dn += 1
                    elif mpra_row.Value > 0:
                        dd += 1
                    else:
                        du += 1
            except Exception as err:
                print("ERR:", homo_seq[i], parent_seq[i], child_seq[i], file=stderr)
                print("ERR:", mpra_data.loc[args.enhancer_name, homo_pos], file=stderr)
                print("ERR:", err, file=stderr)

        overall_du += du
        overall_dd += dd
        overall_dn += dn
        print(node.incident_edges()[-1].label, parent_name, node_name, du, dn, "<-- UP")
        print(
            node.incident_edges()[-1].label, parent_name, node_name, dd, dn, "<-- DOWN"
        )

    pu = len(mpra_data.query("Value > 0 and pval < .05"))
    pd = len(mpra_data.query("Value < 0 and pval < .05"))
    pn = len(mpra_data.query("pval > .05"))

    kukn_fisher = fisher_exact([[overall_du, overall_dn], [pu, pn]])
    kdkn_fisher = fisher_exact([[overall_dd, overall_dn], [pd, pn]])

    if pu > 0 and pn > 0 and overall_dn > 0 and overall_du > 0:
        print("Overall Ku/Kn", (overall_du / pu) / (overall_dn / pn), kukn_fisher[1])
    else:
        print(
            f"Overall Ku/Kn not well defined: Ku = {overall_du}/{pu}, Kn = {overall_dn}/{pn}"
        )
    if pd > 0 and pn > 0 and overall_dn > 0 and overall_dd > 0:
        print("Overall Kd/Kn", (overall_dd / pd) / (overall_dn / pn), kdkn_fisher[1])
    else:
        print(
            f"Overall Kd/Kn not well defined: Kd = {overall_dd}/{pd}, Kn = {overall_dn}/{pn}"
        )
