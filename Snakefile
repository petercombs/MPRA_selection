include_taxids = "9443 9989"
exclude_taxids = "9606"
# 9443 = Primates
# 9889 = Rodents
# 9606 = Humans

rule blast_sequence:
    input:
        seq="{enhancer}/sequence.fasta",
        taxidnodes="Reference/nodes.dmp",
    output:
        "{enhancer}/localblast_withdups.fasta"
    shell: """
    module load blast
    blastn -db refseq_genomic \
        -query {input.seq} \
        -outfmt "6 sscinames saccver staxids pident sseq" \
      | python FilterTaxID.py \
            --taxid-column 3 \
            --include {include_taxids} \
            --exclude {exclude_taxids} \
            --output-fasta \
            {input.taxidnodes} - \
      | cat - <(perl -pe 's/>.*/>Homo_sapiens original sequence/' {input.seq})
      > {output}

    """

rule dedup_blast:
    input:
        seqs="{enhancer}/{file}_withdups.fasta"
    output:
        seqs="{enhancer}/{file}.fasta"
    run:
        from Bio import SeqIO
        default_score = 0
        best_recs = {}
        for rec in SeqIO.parse(input.seqs, 'fasta'):
            score = float(rec.description.split()[-1])
            rec.description = rec.id
            if best_recs.get(rec.id, (default_score, None))[0] < score:
                best_recs[rec.id] = (score, rec)
        SeqIO.write([r[1] for r in best_recs.values()], output.seqs, 'fasta')



rule clustalo_align:
    input:
        "{enhancer}/localblast.fasta",
    output:
        "{enhancer}/clustalo.fasta",
    shell: "clustalo --force -i {input} -o {output} -v"

rule get_phylogeny:
    input:
        seqs="{enhancer}/localblast.fasta",
        tree="Reference/nature05634-s2-revised.txt",
    output:
        newick="{enhancer}/mammals.tree",
        nexus="{enhancer}/mammals.nexus",
    conda: "envs/conda.env"
    run:
        from Bio import SeqIO
        from dendropy import Tree

        species_replacements = [
            ('Nomascus_leucogenys', 'Hylobates_gabriellae'),
            ('Piliocolobus_tephrosceles', 'Procolobus_verus'),
            ('Rhinopithecus_bieti', 'Pygathrix_bieti'),
            ('Rhinopithecus_roxellana', 'Pygathrix_roxellana'),
            ('Pongo_abelii', 'Pongo_pygmaeus'),
            ('Papio_anubis', 'Papio_hamadryas'),
            ('Propithecus_coquereli', 'Propithecus_tattersalli'),
            ('Carlito_syrichta','Tarsius_syrichta'),
            ('Chlorocebus_sabaeus', 'Cercopithecus_mitis'),
            ('Cercocebus_atys', 'Cercocebus_galeritus'),
            ('Fukomys_damarensis', 'Cryptomys_damarensis'),
        ]

        recs = {rec.id: rec for rec in SeqIO.parse(input.seqs, 'fasta')}
        species = {recid for recid in recs if not recid[-1].isnumeric()}
        for from_spec, to_spec in species_replacements:
            species.remove(from_spec)
            species.add(to_spec)
        tree = Tree.get_from_path(input.tree, schema='nexus')

        tree.retain_taxa_with_labels([s.replace('_', ' ') for s in species])
        tree.scale_edges(1/91.6) # Approximate age when rodents diverged from primates
        for to_spec, from_spec in species_replacements:
            n = tree.find_node_with_taxon_label(from_spec.replace('_', ' '))
            n.taxon.label = to_spec
        print(tree.as_string('newick').replace('*','').replace('"', '').replace("'", ''), file=open(output.newick, 'w'))
        tree.taxon_namespace.clear()
        tree.update_taxon_namespace()
        print(tree.as_string('nexus').replace('*','').replace('"', '').replace("'", ''), file=open(output.nexus, 'w'))


