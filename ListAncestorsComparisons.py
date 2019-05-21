"""ListAncestorsComparisons

Perform comparisons of ancestral and derived sequences and look for enrichments
in upregulating mutations and downregulating mutations, compared to chance.
Here, up- and down-regulating is defined as the effect that mutation has, on its
own, in the MPRA data provided.
"""

from os import path
from sys import stderr, stdout
from argparse import ArgumentParser, ArgumentTypeError
from numpy.random import permutation
import numpy as np
import pandas as pd
from matplotlib.pyplot import hist, subplot, title, savefig, vlines, tight_layout, ylim
from bisect import bisect
from Bio import SeqIO, AlignIO
from scipy.stats import fisher_exact
from tqdm import tqdm
from dendropy import Tree, Taxon
from dendropy.utility.error import SeedNodeDeletionException

TREE_LABELS = ["SKIP"] + [
    c * i for i in range(1, 10) for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
]


def parse_args():
    "Parse args from command line"

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


def get_mpra_data(args):
    """Load in MPRA data and convert to standard format


    Returns: a data frame with a MultiIndex with [Element, Position, AltBase],
    and columns Value and pval.

        where Position is a 1-based coordinate starting from the start of the
        sequence (similar to the Patwardhan coordinates, rather than absolute
        coordinates on the chromosome)
    """

    in_data = pd.read_csv(args.mpra_data, sep="\t", header=args.header)
    in_data = in_data[in_data.iloc[:, args.element_column] == args.enhancer_name]
    in_data["pos"] = (
        in_data.iloc[:, args.position_column]
        - in_data.iloc[:, args.position_column].min()
        + 1
    )
    in_data.iloc[:, args.alt_column].name = 'Alt'
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
    return mpra_data


def relabel_tree(tree, target_species_fasta):
    """Update tree with branch labels according to the comparisons we expect to
    do

    """
    if target_species_fasta:
        target_species = {
            rec.id.replace("_", " ")
            for rec in SeqIO.parse(target_species_fasta, "fasta")
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
            else:
                pass
            outgroup_root.taxon = tree.taxon_namespace.new_taxon("Outgroup")
            outgroup_root.label = "Outgroup"
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
        edge.label = TREE_LABELS[i]
        i += 1

    for node in tree:
        node.annotations["comparison"] = node.incident_edges()[-1].label
        if target_species_fasta and node is outgroup_root:
            node.annotations["species"] = "Outgroup"
        else:
            node.annotations["species"] = node.taxon.label if node.taxon else node.label


def score_tree(
    tree, alignment, alignment_posns, mpra_data, enhancer_name, verbose=True
):
    """ Score Ku/Kn and Kd/Kn for the given tree

    """

    overall_du = 0
    overall_dd = 0
    overall_dn = 0

    for node in tree:
        if node.parent_node is None:
            continue
        parent_name = node.parent_node.label.replace(" ", "_")
        node_name = (node.label or node.taxon.label).replace(" ", "_")
        if node_name == "Outgroup":
            continue

        align_len = len(alignment[0])
        homo_pos = 0
        homo_seq = alignment[alignment_posns["original"]]
        parent_seq = alignment[alignment_posns[parent_name]]
        child_seq = alignment[alignment_posns[node_name]]
        branch_du = 0
        branch_dd = 0
        branch_dn = 0
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
                    if len(mpra_data.loc[enhancer_name, homo_pos]) < 3:
                        continue
                    mpra_row = mpra_data.loc[enhancer_name, homo_pos, child_seq[i]]
                    if mpra_row.pval > .05:
                        branch_dn += 1
                    elif mpra_row.Value < 0:
                        branch_dd += 1
                    else:
                        branch_du += 1
                elif child_seq[i] == homo_seq[i] and parent_seq[i] != child_seq[i]:
                    if len(mpra_data.loc[enhancer_name, homo_pos]) < 3:
                        raise ValueError(f"Missing one or more bases at {homo_pos}")
                    mpra_row = mpra_data.loc[enhancer_name, homo_pos, parent_seq[i]]
                    if mpra_row.pval > .05:
                        branch_dn += 1
                    elif mpra_row.Value > 0:
                        branch_dd += 1
                    else:
                        branch_du += 1
            except Exception as err:
                print("ERR:", homo_seq[i], parent_seq[i], child_seq[i], file=stderr)
                print("ERR:", mpra_data.loc[enhancer_name, homo_pos], file=stderr)
                print("ERR:", err, file=stderr)

        node.annotations["udn"] = f"{branch_du}U {branch_dn}N {branch_dd}D"
        overall_du += branch_du
        overall_dd += branch_dd
        overall_dn += branch_dn
        print(
            node.incident_edges()[-1].label,
            parent_name,
            node_name,
            f"{branch_du}U",
            f"{branch_dn}N",
            f"{branch_dd}D",
            sep="\t",
            file=stdout if verbose else open("/dev/null", "w"),
        )

    possible_u = len(mpra_data.query("Value > 0 and pval < .05"))
    possible_d = len(mpra_data.query("Value < 0 and pval < .05"))
    possible_n = len(mpra_data.query("pval > .05"))

    kukn_fisher = fisher_exact([[overall_du, overall_dn], [possible_u, possible_n]])
    kdkn_fisher = fisher_exact([[overall_dd, overall_dn], [possible_d, possible_n]])

    print(
        f"Possible up: {possible_u}", file=stdout if verbose else open("/dev/null", "w")
    )
    print(
        f"Possible neutral: {possible_n}",
        file=stdout if verbose else open("/dev/null", "w"),
    )
    print(
        f"Possible down: {possible_d}",
        file=stdout if verbose else open("/dev/null", "w"),
    )

    if possible_u > 0 and possible_n > 0 and overall_dn > 0 and overall_du > 0:
        print(
            "Overall Ku/Kn {} (p={}; ({}/{})/({}/{}))".format(
                (overall_du / possible_u) / (overall_dn / possible_n),
                kukn_fisher[1],
                overall_du,
                possible_u,
                overall_dn,
                possible_n,
            ),
            file=stdout if verbose else open("/dev/null", "w"),
        )
    else:
        print(
            f"Overall Ku/Kn not well defined: "
            + f"Ku = {overall_du}/{possible_u}, Kn = {overall_dn}/{possible_n}",
            file=stdout if verbose else open("/dev/null", "w"),
        )
    if possible_d > 0 and possible_n > 0 and overall_dn > 0 and overall_dd > 0:
        print(
            "Overall Kd/Kn {} (p={})".format(
                (overall_dd / possible_d) / (overall_dn / possible_n), kdkn_fisher[1]
            ),
            file=stdout if verbose else open("/dev/null", "w"),
        )
    else:
        print(
            f"Overall Kd/Kn not well defined: "
            + f"Kd = {overall_dd}/{possible_d}, Kn = {overall_dn}/{possible_n}",
            file=stdout if verbose else open("/dev/null", "w"),
        )
    return kukn_fisher, kdkn_fisher


def shuffle_mpra(mpra_data):
    """Shuffle data by base that we're switching to

    """
    out = mpra_data.copy()
    for base in out.index.levels[2]:
        out.loc[(slice(None), slice(None), base), :] = permutation(
            out.loc[(slice(None), slice(None), base), :]
        )

    return out


if __name__ == "__main__":
    ARGS = parse_args()
    TREE = Tree.get_from_path(ARGS.input_tree, "newick")

    ALIGNMENT = AlignIO.read(ARGS.seqs, "fasta")
    AlignIO.write(
        ALIGNMENT, path.join(path.dirname(ARGS.seqs), "aln.clustal"), "clustal"
    )
    ALIGNMENT_POSNS = {rec.id: i for i, rec in enumerate(ALIGNMENT)}

    MPRA_DATA = get_mpra_data(ARGS)

    relabel_tree(TREE, ARGS.target_species_fasta)

    REAL_DATA = score_tree(
        TREE, ALIGNMENT, ALIGNMENT_POSNS, MPRA_DATA, ARGS.enhancer_name
    )

    if ARGS.output_tree:
        TREE.write_to_path(ARGS.output_tree, "nexus", suppress_annotations=False)

    SHUFFLED_KUKNS = []
    SHUFFLED_KDKNS = []
    SHUFFLED_KUKN_PS = []
    SHUFFLED_KDKN_PS = []
    for i in tqdm(range(1000)):
        SHUFFLED_MPRA = shuffle_mpra(MPRA_DATA)
        KUKN_SHUF, KDKN_SHUF = score_tree(
            TREE,
            ALIGNMENT,
            ALIGNMENT_POSNS,
            SHUFFLED_MPRA,
            ARGS.enhancer_name,
            verbose=False,
        )
        SHUFFLED_KUKNS.append(KUKN_SHUF[0])
        SHUFFLED_KDKNS.append(KDKN_SHUF[0])
        SHUFFLED_KUKN_PS.append(KUKN_SHUF[1])
        SHUFFLED_KDKN_PS.append(KDKN_SHUF[1])
    SHUFFLED_KUKNS.sort()
    SHUFFLED_KDKNS.sort()

    print(
        "Empirical KuKn p-value",
        bisect(SHUFFLED_KUKNS, REAL_DATA[0][0]) / len(SHUFFLED_KUKNS),
    )
    print(
        "Empirical KdKn p-value",
        bisect(SHUFFLED_KDKNS, REAL_DATA[1][0]) / len(SHUFFLED_KDKNS),
    )

    subplot(2, 2, 1)
    hist(np.log2(SHUFFLED_KDKNS), bins=50)
    vlines(np.log2(REAL_DATA[1][0]), *ylim())
    title("log2 Kd/Kn shuffled")
    subplot(2, 2, 2)
    hist(np.log2(SHUFFLED_KUKNS), bins=50)
    vlines(np.log2(REAL_DATA[0][0]), *ylim())
    title("log2 Ku/Kn shuffled")
    subplot(2, 2, 3)
    hist(SHUFFLED_KDKN_PS, bins=50)
    title("Kd/Kn shuffled p-values")
    subplot(2, 2, 4)
    hist(SHUFFLED_KUKN_PS, bins=50)
    title("Ku/Kn shuffled p-values")

    tight_layout()

    savefig(ARGS.output_tree + ".shuffledpvals.png")
