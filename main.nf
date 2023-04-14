nextflow.enable.dsl=2

import groovyx.gpars.dataflow.DataflowVariable

params.input_fasta = "$baseDir/data/sequences.fasta"
params.outdir = "nf_results"

params.group_by = "country year month"
params.sequences_per_group = 20
params.min_date = 2012
params.dropped_strains = "config/dropped_strains.txt"
params.input_metadata = "data/metadata.tsv"

params.reference = "config/zika_outgroup.gb"

params.coalescent = "opt"
params.date_inference = "marginal"
params.clock_filter_iqd = 4

params.inference_method = "joint"


process index_sequences {
    tag "index_sequences"
    publishDir "${params.outdir}/sequence_index", mode: 'copy'

    input:
    path sequences

    output:
    path "sequence_index.tsv"

    script:
    """
    augur index \
        --sequences $sequences \
        --output sequence_index.tsv
    """
}

process filter_sequences {
    tag "filter_sequences"
    publishDir "${params.outdir}/filtered_sequences", mode: 'copy'

    input:
    path sequences
    path sequence_index
    path metadata
    path exclude

    output:
    path "filtered.fasta"

    script:
    """
    augur filter \
        --sequences $sequences \
        --sequence-index $sequence_index \
        --metadata $metadata \
        --exclude $exclude \
        --output filtered.fasta \
        --group-by ${params.group_by} \
        --sequences-per-group ${params.sequences_per_group} \
        --min-date ${params.min_date}
    """
}

process align_sequences {
    tag "align_sequences"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    path sequences
    path reference

    output:
    path "aligned.fasta"

    script:
    """
    augur align \
        --sequences $sequences \
        --reference-sequence $reference \
        --output aligned.fasta \
        --fill-gaps
    """
}

process build_tree {
    tag "build_tree"
    publishDir "${params.outdir}/tree", mode: 'copy'

    input:
    path alignment
    
    output:
    path "tree_raw.nwk"

    script:
    """
    augur tree \
        --alignment $alignment \
        --output tree_raw.nwk
    """
}

process refine_tree {
    tag "refine_tree"
    publishDir "${params.outdir}/refined_tree", mode: 'copy'

    input:
    path tree
    path alignment
    path metadata
    val coalescent
    val date_inference
    val clock_filter_iqd

    output:
    path "tree_refined.nwk"
    path "branch_lengths.json"

    script:
    """
    augur refine \
        --tree $tree \
        --alignment $alignment \
        --metadata $metadata \
        --output-tree tree_refined.nwk \
        --output-node-data branch_lengths.json \
        --timetree \
        --coalescent $coalescent \
        --date-confidence \
        --date-inference $date_inference \
        --clock-filter-iqd $clock_filter_iqd
    """
}

process ancestral_sequences {
    tag "ancestral_sequences"
    publishDir "${params.outdir}/ancestral_sequences", mode: 'copy'

    input:
    tuple val(method), path(tree), path(alignment)

    output:
    path "nt_muts.json"

    script:
    """
    augur ancestral \
        --tree $tree \
        --alignment $alignment \
        --output-node-data nt_muts.json \
        --inference $method
    """
}


Channel.fromPath(params.input_fasta, checkIfExists:true)
    .set{ input_fasta_ch }
Channel.fromPath(params.input_metadata, checkIfExists:true)
    .set { input_metadata_ch }
Channel.fromPath(params.dropped_strains, checkIfExists:true)
    .set { dropped_strains_ch }
Channel.fromPath(params.reference, checkIfExists:true)
    .set { reference_ch }
// sixth process
Channel
    .value(params.inference_method)
    .set{ inference_method_ch }


/*workflow {
    index_sequences(input_fasta_ch)
    filter_sequences(input_fasta_ch, index_sequences.out, input_metadata_ch, dropped_strains_ch)
    align_sequences(filter_sequences.out, reference_ch)
    align_sequences_ch(filtered_sequences_ch)
        .set{ alignment_ch }
    build_tree(alignment_ch)
}*/

workflow {
    index_sequences(input_fasta_ch)
        .set{ sequence_index_ch }

    filter_sequences(input_fasta_ch, sequence_index_ch, input_metadata_ch, dropped_strains_ch)
        .set{ filtered_sequences_ch }

    align_sequences(filtered_sequences_ch, reference_ch)
        .set{ alignment_ch }

    build_tree(alignment_ch)
        .set{ tree_ch }

    refine_tree(tree_ch, alignment_ch, input_metadata_ch, params.coalescent, params.date_inference, params.clock_filter_iqd)
        .set{ refined_tree_ch }

    ancestral_sequences_ch = Channel
        .from(refined_tree_ch)
        .cross(alignment_ch)
        .map { tree, alignment -> tuple(tree, alignment) }
        .map { items -> tuple(params.inference_method, items[0], items[1]) }

    ancestral_sequences(ancestral_sequences_ch)
}

