from argparse import ArgumentParser, FileType
from tqdm import tqdm
from sys import stderr, stdin
import textwrap as tw


class TaxonomicNode:
    """A taxonomic node

    Various utility functions
    """

    def __init__(self, taxid, parent):
        self.parent = parent
        self.children = []
        self.taxid = taxid

    def is_descendent(self, taxids):
        """Return true if this node is directly beneath any of the given
        taxids"""

        if self.parent is None or not taxids:
            # We are at the root of the tree (or we were provided an empty list)
            return False
        elif self.taxid in taxids:
            # This node is one of the targets
            return True
        else:
            # This node isn't one of the targets, but maybe its parent is?
            return self.parent.is_descendent(taxids)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--include",
        "-i",
        nargs="*",
        default=set(),
        help="""Require each entry to be a child of one (or more) of these taxids""",
    )
    parser.add_argument(
        "--exclude",
        "-x",
        nargs="*",
        default=set(),
        help="""Exclude entries that are children of any of these taxids""",
    )
    parser.add_argument(
        "--taxid-column",
        "-k",
        type=int,
        default=None,
        help="1-based index of the column with the taxids",
    )
    parser.add_argument("--delimiter", "-d", default="\t")
    parser.add_argument("--output-fasta", "-F", default=False, action="store_true")
    parser.add_argument("nodes", help="The nodes.dmp file from NCBI")
    parser.add_argument("input", type=FileType("r"))

    args = parser.parse_args()

    if args.taxid_column is not None:
        args.taxid_column -= 1

    if (not args.include) and (not args.exclude):
        raise ValueError("Must provide at least one include or exclude entry")

    return args


def parse_nodefile(nodes):
    all_nodes = {}
    print("Loading file", file=stderr)
    for line in tqdm(open(nodes)):
        data = [d.strip() for d in line.strip().split("|")]
        taxid = data[0]
        parentid = data[1]
        if taxid == parentid:
            parent = None
        else:
            parent = all_nodes.setdefault(parentid, TaxonomicNode(parentid, None))

        entry = all_nodes.setdefault(taxid, TaxonomicNode(taxid, parent))

        # If this entry was created as the parent of another node then it won't
        # have a parent loaded in by default
        entry.parent = parent

    num_no_parents = 0
    for node in all_nodes.values():
        num_no_parents += node.parent is None
    print(num_no_parents, file=stderr)
    return all_nodes


if __name__ == "__main__":
    args = parse_args()
    nodes = parse_nodefile(args.nodes)

    #print("Hsap in primates", nodes["9606"].is_descendent(["9443"]))
    #print("Horse in primates", nodes["9796"].is_descendent(["9443"]))

    for line in args.input:
        if args.taxid_column is not None:
            columns = line.rstrip("\n").split(args.delimiter)
            data = columns[args.taxid_column].strip()
        else:
            data = line.rstrip("\n")

        if data not in nodes:
            print(
                "Filtered due to not being in the database: ", line.strip(), file=stderr
            )
            continue
        node = nodes[data]
        if (not node.is_descendent(args.exclude)) and node.is_descendent(args.include):
            if args.output_fasta and args.taxid_column is not None:
                species = "_".join(columns[0].split()[:2])
                print(
                    ">{} {}\n{}".format(
                        species,
                        " ".join(columns[1:-1]),  # Other data
                        "\n".join(tw.wrap(columns[-1].replace('-', ''), 100, break_on_hyphens=False)),
                        # Sequence wrapped to 100
                    )
                )
            else:
                print(line, end="")
        elif not node.is_descendent(args.exclude):
            print("Filtered due to not being included: ", line.strip(), file=stderr)
        elif node.is_descendent(args.include):
            print("Filtered due to being excluded: ", line.strip(), file=stderr)
        else:
            print("Should not make it here: ", line.strip(), file=stderr)
