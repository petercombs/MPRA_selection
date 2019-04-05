primates = "9443"
rodents = "9889"
human = "9606"
include_taxids = "{primates} {rodents}"
exclude_taxids = "{human}"

rule reduce_to_mammals:
    input: "Reference/nodes.dmp"
    output: "Reference/mammals.dmp"
    shell: """
    python FilterTaxID.py \
        --delimiter "|" \
        --include 40674 \
        -k 1 \
        {input} {input} \
        > {output}
    """

rule blast_sequence:
    input:
        seq="{enhancer}/sequence.fasta",
    output:
        "{enhancer}/alllocalblast.tsv"
    shell: """
    module load blast
    blastn -db refseq_genomic \
        -query {input.seq} \
        -outfmt "6 sscinames saccver staxids pident sseq" \
      > {output}

    """

rule filter_blast:
    input:
        blastout = "{enhancer}/alllocalblast.tsv",
        seq="{enhancer}/sequence.fasta",
        taxidnodes="Reference/mammals.dmp",
    output:
        "{enhancer}/localblast_withdups.fasta"
    shell:"""
      cat {input.blastout} \
      | python FilterTaxID.py \
            --taxid-column 3 \
            --include {include_taxids} \
            --exclude {exclude_taxids} \
            --output-fasta \
            {input.taxidnodes} - \
      | cat - <(perl -pe 's/>.*/>Homo_sapiens original sequence 100/' {input.seq}) \
      > {output}
      """

rule primates_fasta:
    input:
        blastout = "{enhancer}/alllocalblast.tsv",
        seq="{enhancer}/sequence.fasta",
        taxidnodes="Reference/mammals.dmp",
    output:
        "{enhancer}/primates_withdups.fasta"
    shell:"""
      cat {input.blastout} \
      | python FilterTaxID.py \
            --taxid-column 3 \
            --include {primates} \
            --output-fasta \
            {input.taxidnodes} - \
      | cat - <(perl -pe 's/>.*/>Homo_sapiens original sequence 100/' {input.seq}) \
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
