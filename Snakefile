primates = "9443"
rodents = "314147" # and rabbits
human = "9606"
include_taxids = "{primates} {rodents}".format(primates=primates, rodents=rodents)
exclude_taxids = "{human}".format(human=human)

configfile: "config.yaml"

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

rule clustalw_align:
    input:
        "{enhancer}/localblast.fasta",
    output:
        "{enhancer}/clustalw.clustal",
    shell: "clustalw -infile={input} -align -outfile={output} "

rule clustal_to_fasta:
    input:
        "{file}.clustal"
    output:
        "{file}.fasta"
    #conda: "envs/conda.env"
    run:
        from Bio import SeqIO
        in_recs = SeqIO.parse(input[0], 'clustal')
        SeqIO.write(in_recs, output[0], 'fasta')

rule muscle_align:
    input:
        "{enhancer}/localblast.fasta",
    output:
        "{enhancer}/muscle.fasta",
    shell: "muscle -in {input} -out {output} -diags "

rule tcoffee_align:
    input:
        "{enhancer}/localblast.fasta",
    output:
        aln="{enhancer}/tcoffee.clustal",
        tree="{enhancer}/tcoffee.tree"
    shell: """
    t_coffee -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

rule mcoffee_align:
    input:
        "{enhancer}/localblast.fasta",
    output:
        aln="{enhancer}/mcoffee.clustal",
        tree="{enhancer}/mcoffee.tree"
    shell: """
    t_coffee \
        -method ktup_msa clustalo_msa clustalw2_msa mafftdef_msa dialigntx_msa muscle_msa t_coffee_msa \
        -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

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

rule fastml_reconstruction:
    input:
        tree="{enhancer}/mammals.leaves.tree",
        seqs="{enhancer}/{aligner}.fasta",
    output:
        isdone="{enhancer}/fastml-{aligner}_done",
        log="{enhancer}/FastML-{aligner}/fastml.std",
        seqs="{enhancer}/FastML-{aligner}/seq.marginal_IndelAndChars.txt",
        tree="{enhancer}/FastML-{aligner}/tree.newick.txt",
    shell:"""
    rm -rf `dirname {output.log}`
    perl tools/FastML.v3.11/www/fastml/FastML_Wrapper.pl \
        --MSA_File $PWD/{input.seqs} \
        --seqType NUC \
        --outdir $PWD/`dirname {output.log}` \
        --Tree $PWD/{input.tree} \
        --indelCutOff 0.9 \
        --optimizeBL no
    touch {output.isdone}
    """

rule sequence_from_file:
    input:
        data=lambda wildcards: config['data_files'][wildcards.enhancer],
        exists="enhancers/{enhancer}/exists",
    output:
        "enhancers/{enhancer}/sequence.fasta"
    shell: """
    python ExtractSequenceFromDataFile.py {wildcards.enhancer} {input.data} > {output}
    """

rule exists:
    output: touch("{dir}/exists")

localrules: exists

rule merge_reconstructions:
    input:
        seqs=expand("{{enhancer}}/FastML-{aligner}/seq.marginal_IndelAndChars.txt",
                aligner=['clustalo', 'clustalw', 'tcoffee', 'mcoffee', 'muscle']),
        trees=expand("{{enhancer}}/FastML-{aligner}/tree.newick.txt",
                aligner=['clustalo', 'clustalw', 'tcoffee', 'mcoffee', 'muscle']),

        outdir_exists="{enhancer}/merged/exists",

    output:
        seq="{enhancer}/merged/merged_seq.fasta",
        tree="{enhancer}/merged/tree.newick.txt",
    conda: "envs/conda.env"
    shell: """
    cp {input.trees[0]} {output.tree}
    python CallConsensus.py {output.seq} {input.seqs}

    """

rule ancestor_comparisons:
    input:
        primates="enhancers/{enhancer}/primates.fasta",
        tree="enhancers/{enhancer}/FastML/tree.newick.txt",
        seq="enhancers/{enhancer}/FastML/seq.marginal_IndelAndChars.txt",
        data=lambda wildcards: config['data_files'][wildcards.enhancer],
        script="ListAncestorsComparisons.py",
    output:
        results="enhancers/{enhancer}/FastML/selection_results.txt",
        tree="enhancers/{enhancer}/FastML/comparisons.tree",
    params:
        ename=lambda wildcards: (wildcards.enhancer.lower()
                if 'Patwardhan' in config['data_files'][wildcards.enhancer]
                else wildcards.enhancer)
    conda: "envs/conda.env"
    shell: """
    python {input.script} \
        --enhancer-name {params.ename} \
        -t {input.primates} \
        --output-tree {output.tree} \
        {input.seq} {input.tree} {input.data} \
        > {output.results}
    """

rule all_selection:
    input:
        expand("enhancers/{enhancer}/{recon}/selection_results.txt",
                enhancer=config["data_files"].keys(),
                recon=["FastML"],
        )
