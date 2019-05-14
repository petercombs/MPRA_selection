mammals = "40674"
primates = "9443"
rodents = "314147" # and rabbits
human = "9606"

groups = {
    "mammals" : "40674",
    "primates" : "9443",
    "rodents" : "314147",
    "human" : "9606",
}

include_taxids = "{primates} {rodents}".format(primates=primates, rodents=rodents)
exclude_taxids = "{human}".format(human=human)

configfile: "config.yaml"

target_fasta_file='localblast'

rule get_nodes:
    output:
        "Reference/nodes.dmp",
        "Reference/names.dmp",
    shell: """
    cd Reference
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz
    """

rule reduce_to_mammals:
    input: "Reference/nodes.dmp"
    output: "Reference/mammals.dmp"
    conda: "envs/conda.yaml"
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
        stripped='{enhancer}/{file}.fasta.stripped',
        seq='{enhancer}/{file}.fasta.masked',
        joined='{enhancer}/{file}_withmasked.fasta'
    conda: "envs/conda.yaml"
    shell: """
    body () {{
        IFS= read -r header
            printf '%s\n' "$header"
            "$@"
    }}

    perl -pe 's/-//g' {input} > {output.stripped}

    RepeatMasker \
        -lib Reference/RepBase24.03.fasta/all.fasta \
        -s {output.stripped}

    if [ `grep -c "no repetitive sequences detected" {output.stripped}.out` -gt 0 ]
    then
    cp {input} {output.seq}
    else
    mv {output.stripped}.masked {output.seq}
    fi

    cat  {output.seq} \
        | sed '/^[^>]/s/N//g' \
        | sed '/^$/d' \
        | cat {input} -  \
        > {output.joined}
    """

rule propagate_masks:
    input:
        aligned="{file}.fasta",
        masked="{file}.fasta.masked",
    output:
        "{file}_maskprop.fasta",
    shell: """
    python PropagateMasks.py {input} {output}
    """

rule blast_sequence:
    input:
        seq="{enhancer}/sequence_withmasked.fasta",
    output:
        "{enhancer}/alllocalblast.tsv"
    conda: "envs/conda.yaml"
    shell: """
    blastn -db refseq_genomic \
        -query {input.seq} \
        -outfmt "6 sscinames saccver staxids pident evalue sseq" \
      > {output}

    """

def get_species(wildcards):
    return config['species'].get(wildcards.enhancer, 'Homo_sapiens')


rule filter_blast:
    input:
        blastout = "{enhancer}/alllocalblast.tsv",
        seq="{enhancer}/sequence.fasta",
        taxidnodes="Reference/mammals.dmp",
    output:
        "{enhancer}/localblast_withdups.fasta"
    params:
        species=get_species,
    conda: "envs/conda.yaml"
    shell:"""
      cat {input.blastout} \
      | python FilterTaxID.py \
            --taxid-column 3 \
            --include {include_taxids} \
            --exclude {exclude_taxids} \
            --output-fasta \
            {input.taxidnodes} - \
      | cat - <(perl -pe 's/>.*/>{params.species} original sequence 100/' {input.seq}) \
      > {output}
      """

rule mammals_fasta:
    input:
        blastout = "enhancers/{enhancer}/alllocalblast.tsv",
        seq="enhancers/{enhancer}/sequence.fasta",
        taxidnodes="Reference/mammals.dmp",
    output:
        "enhancers/{enhancer}/mammals_withdups.fasta"
    params:
        species=get_species,
    conda: "envs/conda.yaml"
    shell:"""
      cat {input.blastout} \
      | python FilterTaxID.py \
            --taxid-column 3 \
            --include {mammals} \
            --output-fasta \
            {input.taxidnodes} - \
      | cat - <(perl -pe 's/>.*/>{params.species} original sequence 100/' {input.seq}) \
      > {output}
      """



rule primates_fasta:
    input:
        blastout = "{enhancer}/alllocalblast.tsv",
        seq="{enhancer}/sequence.fasta",
        taxidnodes="Reference/mammals.dmp",
    output:
        "{enhancer}/{ingroup}_withdups.fasta"
    params:
        ingroup=lambda wildcards: groups[wildcards.ingroup],
        species=get_species,
    conda: "envs/conda.yaml"
    shell:"""
      cat {input.blastout} \
      | python FilterTaxID.py \
            --taxid-column 3 \
            --include {params.ingroup} \
            --output-fasta \
            {input.taxidnodes} - \
      | cat - <(perl -pe 's/>.*/>{params.species} original sequence 100/' {input.seq}) \
      > {output}
      """



rule dedup_blast:
    input:
        seqs="{enhancer}/{file}_withdups.fasta"
    output:
        seqs="{enhancer}/{file}.fasta"
    wildcard_constraints:
        file='[^_]*'
    conda: "envs/conda.yaml"
    shell: "python DedupBlast.py {input} {output}"



rule clustalo_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        "{enhancer}/{target}_clustalo_aligned.fasta",
    conda: "envs/conda.yaml"
    shell: "clustalo --force -i {input} -o {output} -v"

rule clustalw_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        "{enhancer}/{target}_clustalw_aligned.clustal",
    conda: "envs/conda.yaml"
    shell: "clustalw -infile={input} -align -outfile={output} "

rule clustal_to_fasta:
    input:
        "{file}_aligned.clustal"
    output:
        "{file}_aligned.fasta"
    conda: "envs/conda.yaml"
    shell: "python ClustalToFasta.py {input} {output}"

rule muscle_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        "{enhancer}/{target}_muscle_aligned.fasta",
    conda: "envs/conda.yaml"
    shell: "muscle -in {input} -out {output} -diags "

rule tcoffee_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        aln="{enhancer}/{target}_tcoffee_aligned.clustal",
        tree="{enhancer}/{target}_tcoffee.tree"
    conda: "envs/conda.yaml"
    shell: """
    t_coffee -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

rule merged_tcoffee_align:
    input:
        "{enhancer}/{merged}/merged_seq.fasta",
    output:
        aln="{enhancer}/{merged}/merged_aligned.clustal",
        tree="{enhancer}/{merged}/tcoffee.tree"
    conda: "envs/conda.yaml"
    shell: """
    t_coffee -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

rule mcoffee_align:
    input:
        "{enhancer}/{target}_withsupport.fasta",
    output:
        aln="{enhancer}/{target}_mcoffee_aligned.clustal",
        tree="{enhancer}/{target}_mcoffee.tree"
    conda: "envs/conda.yaml"
    shell: """
    t_coffee \
        -method ktup_msa clustalo_msa clustalw2_msa mafftdef_msa dialigntx_msa muscle_msa t_coffee_msa \
        -n_core 10 \
        -seq {input} -outfile {output.aln} -newtree {output.tree}
    """

rule get_phylogeny:
    input:
        seqs="enhancers/{enhancer}/{target}.fasta",
        tree="Reference/nature05634-s2-revised.txt",
        primates=lambda wildcards: "enhancers/{enhancer}/{ingroup}.fasta".format(
                ingroup=config['ingroup'].get(wildcards.enhancer, 'primates'),
                enhancer=wildcards.enhancer),
    output:
        newickunscaled="enhancers/{enhancer}/{target}.unscaled.tree",
        newick="enhancers/{enhancer}/{target}.tree",
        nexus="enhancers/{enhancer}/{target}.nexus",
        fasta="enhancers/{enhancer}/{target}_withsupport.fasta",
    conda: "envs/conda.yaml"
    shell: """
    python GetPhylogeny.py \
        --targets {input.primates} \
        {input.seqs} {input.tree} enhancers/{wildcards.enhancer}/{wildcards.target}
    """

rule strip_internal_nodes:
    input:
        "{tree}.tree"
    output:
        "{tree}.leaves.tree"
    conda: "envs/conda.yaml"
    shell: """
    sed 's/[0-9][0-9_A-Za-z]*:/:/g' < {input} \
        | sed 's/\[&R\] //' \
        > {output}
    """

rule fastml_reconstruction:
    input:
        tree="{enhancer}/{target}.leaves.tree",
        seqs="{enhancer}/{target}_{aligner}_aligned.fasta",
    output:
        isdone="{enhancer}/fastml-{target}-{aligner}_done",
        log="{enhancer}/FastML-{target}-{aligner}/fastml.std",
        seqs="{enhancer}/FastML-{target}-{aligner}/seq.marginal_IndelAndChars.txt",
        tree="{enhancer}/FastML-{target}-{aligner}/tree.newick.txt",
    wildcard_constraints:
        target='[^-_]*'
    conda: "envs/conda.yaml"
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
    params:
        species=get_species,
    conda: "envs/conda.yaml"
    shell: """
    python ExtractSequenceFromDataFile.py {wildcards.enhancer} {input.data} {params.species} > {output}
    """

rule exists:
    output: touch("{dir}/exists")


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
    wildcard_constraints:
        target='[^_-]*'
    conda: "envs/conda.yaml"
    shell: """
    cp {input.trees[0]} {output.tree}
    python CallConsensus.py {output.seq} {input.seqs}
    echo ">original" >> {output.seq}
    awk "NR > 1" {input.original} >> {output.seq}
    """

#ruleorder: sequence_from_file > primates_fasta

ruleorder: sequence_from_file > merge_reconstructions > dedup_blast > mammals_fasta > primates_fasta > repmask_input > propagate_masks > filter_blast >  get_phylogeny > clustalo_align > clustal_to_fasta > muscle_align

ruleorder: strip_internal_nodes > get_phylogeny

rule ancestor_comparisons:
    input:
        primates=lambda wildcards: "enhancers/{enhancer}/{ingroup}.fasta".format(
            ingroup=config['ingroup'].get(wildcards.enhancer, 'primates'),
            enhancer=wildcards.enhancer,
        ),
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
    conda: "envs/conda.yaml"
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
        primates=lambda wildcards: "enhancers/{enhancer}/{ingroup}.fasta".format(
            ingroup=config['ingroup'].get(wildcards.enhancer, 'primates'),
            enhancer=wildcards.enhancer,
        ),
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

    conda: "envs/conda.yaml"
    shell: """
    python {input.script} \
        --enhancer-name {params.ename} \
        --{params.data_style} \
        -t {input.primates} \
        --output-tree {output.tree} \
        {input.seq} {input.tree} {input.data} \
        > {output.results}
    """

rule comparisons_figure:
    input:
        tree="{dir}/comparisons.tree",
        figtreeconf="Reference/figtree.conf",
    output:
        with_plot="{dir}/comparisons.figtree.tree",
        imageraster="{dir}/comparisons.png",
        imagevector="{dir}/comparisons.pdf",
    shell: """
    cat {input} > {output.with_plot}
    figtree -width 1200 -height 900 -graphic PNG {output.with_plot} {output.imageraster}
    figtree -width 1200 -height 900 -graphic PDF {output.with_plot} {output.imagevector}

    """

rule renamed_figures:
    input:
        "figures/exists",
        fig="enhancers/{enhancer}/mammals/comparisons.{type}",
    output:
        fig="figures/{enhancer}.{type}"
    shell: "cp {input.fig} {output.fig}"

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
        ),
        #expand("enhancers/{enhancer}/smith/selection_results.txt",
                #enhancer=['ECR11', 'ALDOB'],
              #),

rule all_figures:
    input:
        expand("figures/{enhancer}.{types}",
            enhancer=config["data_files"].keys(),
            types=["png", "pdf"],
        )


localrules: exists
localrules: all_mammal_seqs, all_repmasked, all_selection, all_figures, renamed_figures
