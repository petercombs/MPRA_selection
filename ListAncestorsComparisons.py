from argparse import ArgumentParser, ArgumentTypeError
from dendropy import Tree, Taxon
from dendropy.utility.error import SeedNodeDeletionException
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
    data_style = parser.add_argument_group(
        "MPRA Data Format",
        (
            "Either choose from a pre-specified data format or specify all the"
            " important columns"
        ),
    )
    data_style.add_argument(
        "--patwardhan",
        dest="data_style",
        action="store_const",
        const="patwardhan",
        help="Data format matches Patwardhan et al 2012 (Nat Biotech)",
    )
    data_style.add_argument(
        "--kirchner",
        dest="data_style",
        action="store_const",
        const="kirchner",
        help="Data format matches Kirchner et al 2018 (bioRxiv)",
    )
    data_style.add_argument("--alt-column", type=int)
    data_style.add_argument("--header", action="store_const", default=None, const=0)
    data_style.add_argument(
        "--no-header", dest="header", action="store_const", const=None
    )
    data_style.add_argument("--element-column", type=int)
    data_style.add_argument("--position-column", type=int)
    data_style.add_argument("--value-column", type=int)
    data_style.add_argument("--pval-column", type=int)

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

    if args.data_style == "patwardhan":
        args.header = None
        args.element_column = 0
        args.position_column = 1
        args.alt_column = 2
        args.value_column = 3
        args.pval_column = 4
    elif args.data_style == "kirchner":
        args.header = 0
        args.element_column = 9
        args.position_column = 1
        args.alt_column = 3
        args.value_column = 7
        args.pval_column = 8
    elif None in (
        args.element_column,
        args.position_column,
        args.alt_column,
        args.value_column,
        args.pval_column,
    ):
        raise ArgumentTypeError(
            "One of the data columns has not been specified:"
            + str(
                {
                    "Element column": args.element_column,
                    "Position column": args.position_column,
                    "Alternate base column": args.alt_column,
                    "MPRA Value Column": args.value_column,
                    "P-value column": args.pval_column,
                }
            )
        )

    if args.header is False:
        args.header = None
    return args


if __name__ == "__main__":
    args = parse_args()
    tree = Tree.get_from_path(args.input_tree, "newick")

    alignment = AlignIO.read(args.seqs, "fasta")
    AlignIO.write(
        alignment, path.join(path.dirname(args.seqs), "aln.clustal"), "clustal"
    )
    alignment_posns = {rec.id: i for i, rec in enumerate(alignment)}

    """
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
            data={
                "ref": list(mpra_data["Ref"]),
                "OldPos": list(mpra_data["Position"]),
                "Value": list(mpra_data.Value),
                "pval": list(mpra_data["P-Value"]),
            },
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
        """
    in_data = pd.read_csv(args.mpra_data, sep="\t", header=args.header)
    in_data = in_data[in_data.iloc[:, args.element_column] == args.enhancer_name]
    in_data["pos"] = (
        in_data.iloc[:, args.position_column]
        - in_data.iloc[:, args.position_column].min()
        + 1
    )
    index = pd.MultiIndex.from_arrays(
        [
            in_data.iloc[:, args.element_column],
            in_data.pos,
            in_data.iloc[:, args.alt_column],
        ]
    )

    mpra_data = pd.DataFrame(
        index=index,
        data={
            "Value": list(in_data.iloc[:, args.value_column]),
            "pval": list(in_data.iloc[:, args.pval_column]),
        },
    )

    mpra_data.loc[args.enhancer_name].to_csv(
        path.join(path.dirname(args.input_tree), "mpra_data.tsv"), sep="\t"
    )
    print(tree.as_string("newick"))
    if args.target_species_fasta:
        target_species = {
            rec.id.replace("_", " ")
            for rec in SeqIO.parse(args.target_species_fasta, "fasta")
        }
        try:
            outgroup_tree = tree.extract_tree_without_taxa_labels(target_species)
            outgroup_root = tree.find_node(
                lambda n: n.label == outgroup_tree.nodes()[0].label
            )
            if outgroup_root.child_nodes():
                # There's some kind of error that happens if there's only one
                # outgroup species
                outgroup_root.clear_child_nodes()
                outgroup_root.taxon = Taxon(outgroup_root.label)
            else:
                pass
        except SeedNodeDeletionException:
            print(
                "WARNING: No species were provided that are not in the target group",
                file=stderr,
            )
            outgroup_root = None
    else:
        outgroup_root = None

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
            node.annotations["species"] = "Outgroup"
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
        homo_seq = alignment[alignment_posns["original"]]
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
        print(
            node.incident_edges()[-1].label,
            parent_name,
            node_name,
            f"{du}U",
            f"{dn}N",
            f"{dd}D",
            sep="\t",
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
