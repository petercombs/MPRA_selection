from sys import stderr
from Bio import SeqIO
from Bio.Seq import Seq
from dendropy import Tree

species_replacements = [
    ("Nomascus_leucogenys", "Hylobates_gabriellae"),
    ("Piliocolobus_tephrosceles", "Procolobus_verus"),
    ("Rhinopithecus_bieti", "Pygathrix_bieti"),
    ("Rhinopithecus_roxellana", "Pygathrix_roxellana"),
    ("Pongo_abelii", "Pongo_pygmaeus"),
    ("Papio_anubis", "Papio_hamadryas"),
    ("Propithecus_coquereli", "Propithecus_tattersalli"),
    ("Carlito_syrichta", "Tarsius_syrichta"),
    ("Chlorocebus_sabaeus", "Cercopithecus_mitis"),
    ("Cercocebus_atys", "Cercocebus_galeritus"),
    ("Fukomys_damarensis", "Cryptomys_damarensis"),
    ("Ictidomys_tridecemlineatus", "Spermophilus_tridecemlineatus"),
    ("Cricetulus_griseus", "Cricetulus_barabensis"),
    ("Urocitellus_parryii", "Spermophilus_parryii"),
    ("Nannospalax_galili", "Nannospalax_ehrenbergi"),
    ("Plecturocebus_moloch", "Callicebus_moloch"),
]

replacement_dict = {b: a for a, b in species_replacements}


def parse_args():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--targets", default=False)
    parser.add_argument("seqs")
    parser.add_argument("in_tree")
    parser.add_argument("out_base_name")
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = parse_args()

    recs = {rec.id: rec for rec in SeqIO.parse(args.seqs, "fasta")}

    empty_seqs = []
    for recid in recs:
        recs[recid].seq = Seq(str(recs[recid].seq).replace("N", "").replace("-", ""))
        if len(recs[recid].seq) == 0:
            empty_seqs.append(recid)
    for recid in empty_seqs:
        recs.pop(recid)

    species = {recid for recid in recs if not recid[-1].isnumeric()}
    for from_spec, to_spec in species_replacements:
        if from_spec in species:
            species.remove(from_spec)
            species.add(to_spec)
    tree = Tree.get_from_path(args.in_tree, schema="nexus")
    tree_taxa = {tx.label.replace(" ", "_") for tx in tree.taxon_namespace}

    print(species.difference(tree_taxa), file=stderr)

    tree.retain_taxa_with_labels([s.replace("_", " ") for s in species])
    if args.targets:
        orig_tree = Tree(tree)
        targets = {rec.id for rec in SeqIO.parse(args.targets, "fasta")}
        target_labels = [
            l.replace("_", " ")
            for l in targets.intersection(tree_taxa).intersection(species)
        ]
        target_root = tree.mrca(taxon_labels=target_labels)
        tree = Tree(seed_node=(target_root.parent_node or target_root))

    tree_taxa = {tx.label.replace(" ", "_") for tx in tree.taxon_namespace}
    for to_spec, from_spec in species_replacements:
        n = tree.find_node_with_taxon_label(from_spec.replace("_", " "))
        if n:
            n.taxon.label = to_spec
    print(
        tree.as_string("newick").replace("*", "").replace('"', "").replace("'", ""),
        file=open(args.out_base_name + ".unscaled.tree", "w"),
    )
    tree.scale_edges(1 / 91.6)  # Approximate age when rodents diverged from primates
    print(
        tree.as_string("newick").replace("*", "").replace('"', "").replace("'", ""),
        file=open(args.out_base_name + ".tree", "w"),
    )
    tree.taxon_namespace.clear()
    tree.update_taxon_namespace()
    print(
        tree.as_string("nexus").replace("*", "").replace('"', "").replace("'", ""),
        file=open(args.out_base_name + ".nexus", "w"),
    )

    SeqIO.write(
        [
            recs[replacement_dict.get(spec, spec)]
            for spec in species.intersection(tree_taxa)
        ],
        open(args.out_base_name + "_withsupport.fasta", "w"),
        "fasta",
    )
