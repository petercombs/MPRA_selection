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
        newickunscaled="{enhancer}/mammals.unscaled.tree",
        newick="{enhancer}/mammals.tree",
        nexus="{enhancer}/mammals.nexus",
    conda: "envs/conda.env"
    shell: "python GetPhylogeny.py {input} {output}"

rule strip_internal_nodes:
    input:
        "{tree}.tree"
    output:
        "{tree}.leaves.tree"
    shell: """
    sed 's/[0-9][0-9_A-Za-z]*:/:/g' < {input} \
        | sed 's/\[&R\] //' \
        > {output}
    """
