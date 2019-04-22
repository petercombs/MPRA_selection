from sys import stderr, argv
from Bio import SeqIO
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

if __name__ == "__main__":
    recs = {rec.id: rec for rec in SeqIO.parse(argv[1], "fasta")}
    species = {recid for recid in recs if not recid[-1].isnumeric()}
    for from_spec, to_spec in species_replacements:
        if from_spec in species:
            species.remove(from_spec)
            species.add(to_spec)
    tree = Tree.get_from_path(argv[2], schema="nexus")
    tree_taxa = {tx.label.replace(" ", "_") for tx in tree.taxon_namespace}

    print(species.difference(tree_taxa), file=stderr)

    tree.retain_taxa_with_labels([s.replace("_", " ") for s in species])
    for to_spec, from_spec in species_replacements:
        n = tree.find_node_with_taxon_label(from_spec.replace("_", " "))
        if n:
            n.taxon.label = to_spec
    print(
        tree.as_string("newick").replace("*", "").replace('"', "").replace("'", ""),
        file=open(argv[3], "w"),
    )
    tree.scale_edges(1 / 91.6)  # Approximate age when rodents diverged from primates
    print(
        tree.as_string("newick").replace("*", "").replace('"', "").replace("'", ""),
        file=open(argv[4], "w"),
    )
    tree.taxon_namespace.clear()
    tree.update_taxon_namespace()
    print(
        tree.as_string("nexus").replace("*", "").replace('"', "").replace("'", ""),
        file=open(argv[5], "w"),
    )

    SeqIO.write(
        [
            recs[replacement_dict.get(spec, spec)]
            for spec in species.intersection(tree_taxa)
        ],
        open(argv[6], "w"),
        "fasta",
    )
