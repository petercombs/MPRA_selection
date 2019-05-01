mammals = "40674"
primates = "9443"
rodents = "314147" # and rabbits
human = "9606"
include_taxids = "{primates} {rodents}".format(primates=primates, rodents=rodents)
exclude_taxids = "{human}".format(human=human)

configfile: "config.yaml"

target_fasta_file='localblast'
rule reduce_to_mammals:
    input: "Reference/nodes.dmp"
    output: "Reference/mammals.dmp"
    conda: "envs/conda.env"
    shell: """
    python FilterTaxID.py \
        --delimiter "|" \
        --include {mammals} \
        -k 1 \
        {input} {input} \
        > {output}
    """

rule repmask_input:
    input:
        seq='{enhancer}/{file}.fasta'
    output:
        seq='{enhancer}/{file}.fasta.masked',
        joined='{enhancer}/{file}_withmasked.fasta'
    conda: "envs/conda.env"
    shell: """
    body () {{
        IFS= read -r header
            printf '%s\n' "$header"
            "$@"
    }}

    RepeatMasker \
        -lib Reference/RepBase24.03.fasta/all.fasta \
        -s {input}

    if [ `grep -c "no repetitive sequences detected" {input}.out` -gt 0 ]
    then
    cp {input} {output.seq}
    fi

    cat  {output.seq} \
        | body sed 's/N//g' \
        | sed '/^$/d' \
        | cat {input} -  \
        > {output.joined}
    """

rule blast_sequence:
    input:
        seq="{enhancer}/sequence_withmasked.fasta",
    output:
        "{enhancer}/alllocalblast.tsv"
    conda: "envs/conda.env"
    shell: """
    module load blast
    blastn -db refseq_genomic \
        -query {input.seq} \
        -outfmt "6 sscinames saccver staxids pident evalue sseq" \
      > {output}

    """



rule filter_blast:
    input:
        blastout = "{enhancer}/alllocalblast.tsv",
        seq="{enhancer}/sequence.fasta",
        taxidnodes="Reference/mammals.dmp",
    output:
        "{enhancer}/localblast_withdups.fasta"
    conda: "envs/conda.env"
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

rule mammals_fasta:
    input:
        blastout = "{enhancer}/alllocalblast.tsv",
        seq="{enhancer}/sequence.fasta",
        taxidnodes="Reference/mammals.dmp",
    output:
        "{enhancer}/mammals_withdups.fasta"
    conda: "envs/conda.env"
    shell:"""
      cat {input.blastout} \
      | python FilterTaxID.py \
            --taxid-column 3 \
            --include {mammals} \
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
    conda: "envs/conda.env"
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
        default_score = 1000
        best_recs = {}
        for rec in SeqIO.parse(input.seqs, 'fasta'):
            score = float(rec.description.split()[-1])
            rec.description = rec.id
            if best_recs.get(rec.id, (default_score, None))[0] > score:
                best_recs[rec.id] = (score, rec)
        SeqIO.write([r[1] for r in best_recs.values()], output.seqs, 'fasta')



rule clustalo_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        "{enhancer}/{target}_clustalo.fasta",
    conda: "envs/conda.env"
    shell: "clustalo --force -i {input} -o {output} -v"

rule clustalw_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        "{enhancer}/{target}_clustalw.clustal",
    conda: "envs/conda.env"
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
        "{enhancer}/{target}_withsupport.fasta",
    output:
        "{enhancer}/{target}_muscle.fasta",
    conda: "envs/conda.env"
    shell: "muscle -in {input} -out {output} -diags "

rule tcoffee_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        aln="{enhancer}/{target}_tcoffee.clustal",
        tree="{enhancer}/{target}_tcoffee.tree"
    conda: "envs/conda.env"
    shell: """
    t_coffee -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

rule merged_tcoffee_align:
    input:
        "{enhancer}/{merged}/merged_seq.fasta",
    output:
        aln="{enhancer}/{merged}/merged_aligned.clustal",
        tree="{enhancer}/{merged}/tcoffee.tree"
    conda: "envs/conda.env"
    shell: """
    t_coffee -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

rule mcoffee_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        aln="{enhancer}/{target}_mcoffee.clustal",
        tree="{enhancer}/{target}_mcoffee.tree"
    conda: "envs/conda.env"
    shell: """
    t_coffee \
        -method ktup_msa clustalo_msa clustalw2_msa mafftdef_msa dialigntx_msa muscle_msa t_coffee_msa \
        -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

rule get_phylogeny:
    input:
        seqs="{enhancer}/{target}.fasta.masked",
        tree="Reference/nature05634-s2-revised.txt",
        primates="{enhancer}/primates.fasta",
    output:
        newickunscaled="{enhancer}/{target}.unscaled.tree",
        newick="{enhancer}/{target}.tree",
        nexus="{enhancer}/{target}.nexus",
        fasta="{enhancer}/{target}_withsupport.fasta",
    conda: "envs/conda.env"
    shell: """
    python GetPhylogeny.py \
        --targets {input.primates} \
        {input.seqs} {input.tree} {wildcards.enhancer}/{wildcards.target}
    """

rule strip_internal_nodes:
    input:
        "{tree}.tree"
    output:
        "{tree}.leaves.tree"
    conda: "envs/conda.env"
    shell: """
    sed 's/[0-9][0-9_A-Za-z]*:/:/g' < {input} \
        | sed 's/\[&R\] //' \
        > {output}
    """

rule fastml_reconstruction:
    input:
        tree="{enhancer}/{target}.leaves.tree",
        seqs="{enhancer}/{target}_{aligner}.fasta",
    output:
        isdone="{enhancer}/fastml-{target}-{aligner}_done",
        log="{enhancer}/FastML-{target}-{aligner}/fastml.std",
        seqs="{enhancer}/FastML-{target}-{aligner}/seq.marginal_IndelAndChars.txt",
        tree="{enhancer}/FastML-{target}-{aligner}/tree.newick.txt",
    conda: "envs/conda.env"
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
    conda: "envs/conda.env"
    shell: """
    python ExtractSequenceFromDataFile.py {wildcards.enhancer} {input.data} > {output}
    """

rule exists:
    output: touch("{dir}/exists")

localrules: exists

rule merge_reconstructions:
    input:
        original="{enhancer}/sequence.fasta",
        seqs=lambda wildcards: expand("{enhancer}/FastML-{target}-{aligner}/seq.marginal_IndelAndChars.txt",
                enhancer=[wildcards.enhancer],
                target=[wildcards.target],
                aligner=['clustalo', 'clustalw', 'tcoffee', 'mcoffee', 'muscle']),
        trees=lambda wildcards: expand("{enhancer}/FastML-{target}-{aligner}/tree.newick.txt",
                enhancer=[wildcards.enhancer],
                target=[wildcards.target],
                aligner=['clustalo', 'clustalw', 'tcoffee', 'mcoffee', 'muscle']),
        outdir_exists="{enhancer}/{target}/exists",
    output:
        seq="{enhancer}/{target}/merged_seq.fasta",
        tree="{enhancer}/{target}/tree.newick.txt",
    conda: "envs/conda.env"
    shell: """
    cp {input.trees[0]} {output.tree}
    python CallConsensus.py {output.seq} {input.seqs}
    echo ">original" >> {output.seq}
    awk "NR > 1" {input.original} >> {output.seq}
    """

ruleorder: merge_reconstructions > repmask_input > propagate_masks > filter_blast > dedup_blast > clustalo_align > clustal_to_fasta > muscle_align

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

rule merged_ancestor_comparisons:
    input:
        primates="enhancers/{enhancer}/primates.fasta",
        tree="enhancers/{enhancer}/{merged}/tree.newick.txt",
        seq="enhancers/{enhancer}/{merged}/merged_aligned_maskprop.fasta",
        data=lambda wildcards: config['data_files'][wildcards.enhancer],
        script="ListAncestorsComparisons.py",
    output:
        results="enhancers/{enhancer}/{merged}/selection_results.txt",
        tree="enhancers/{enhancer}/{merged}/comparisons.tree",
    params:
        ename=lambda wildcards: (wildcards.enhancer.lower()
                if 'Patwardhan' in config['data_files'][wildcards.enhancer]
                else wildcards.enhancer),
        data_style=lambda wildcards: config['data_styles'][wildcards.enhancer]

    conda: "envs/conda.env"
    shell: """
    python {input.script} \
        --enhancer-name {params.ename} \
        --{params.data_style} \
        -t {input.primates} \
        --output-tree {output.tree} \
        {input.seq} {input.tree} {input.data} \
        > {output.results}
    """

rule all_mammal_seqs:
    input:
        expand("enhancers/{enhancer}/mammals.fasta",
                enhancer=config["data_files"].keys(),
        )

rule all_repmasked:
    input:
        expand("enhancers/{enhancer}/sequence_withmasked.fasta",
                enhancer=config["data_files"].keys(),
        )

rule all_selection:
    input:
        expand("enhancers/{enhancer}/{reconstruction}/selection_results.txt",
                enhancer=config["data_files"].keys(),
                reconstruction=["mammals",],
        )
