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

params.colors = file("config/colors.tsv")
params.lat_longs = file("config/lat_longs.tsv")
params.auspice_config = file("config/auspice_config.json")


params.coalescent = "opt"
params.date_inference = "marginal"
params.clock_filter_iqd = 4

params.inference_method = "joint"
params.columns = "region country"

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
    tuple path("tree_refined.nwk"), path("branch_lengths.json")

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

process translate {
    tag "translate"
    publishDir "${params.outdir}/translated_sequences", mode: 'copy'

    input:
    tuple path(tree), path(node_data), path(reference)

    output:
    path "aa_muts.json"

    script:
    """
    augur translate \
        --tree $tree \
        --ancestral-sequences $node_data \
        --reference-sequence $reference \
        --output-node-data aa_muts.json
    """
}

process traits {
    tag "traits"
    publishDir "${params.outdir}/traits", mode: 'copy'

    input:
        tuple val(columns), path(tree), path(metadata)

    output:
        tuple path("traits.json"), path("traitsregion.mugration_model.txt"), path("traitscountry.mugration_model.txt")

    script:
    """
    augur traits \
        --tree ${tree} \
        --metadata ${metadata} \
        --output-node-data traits.json \
        --columns ${columns} \
        --confidence
    """
}

process export {
    tag "export"
    publishDir "${params.outdir}/export", mode: 'copy'

    input:
        tuple path(tree), path(metadata), path(branch_lengths), path(traits), path(nt_muts), path(aa_muts), path(colors), \
            path(lat_longs), path(auspice_config)

    output:
        path("auspice_json.json")

    log.info "Exporting data files for auspice..."

    script:
        """
        augur export v2 \
            --tree ${tree} \
            --metadata ${metadata} \
            --node-data ${branch_lengths} ${traits} ${nt_muts} ${aa_muts} \
            --colors ${colors} \
            --lat-longs ${lat_longs} \
            --auspice-config ${auspice_config} \
            --include-root-sequence \
            --output auspice_json.json
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
Channel
    .value(params.columns)
    .set{ columns_ch }
Channel
    .value(params.colors)
    .set{ colors_ch }
Channel
    .value(params.lat_longs)
    .set { lat_longs_ch }
Channel 
    .value(params.auspice_config)
    .set { auspice_config_ch }

workflow {
    index_sequences(input_fasta_ch)
        .set{ sequence_index_ch }

    filter_sequences(input_fasta_ch, sequence_index_ch, input_metadata_ch, dropped_strains_ch)
        .set{ filtered_sequences_ch }

    align_sequences(filtered_sequences_ch, reference_ch)
        .set{ alignment_ch }
 //   alignment_ch |view

    build_tree(alignment_ch)
        .set{ tree_ch }

    refine_tree(tree_ch, alignment_ch, input_metadata_ch, params.coalescent, params.date_inference, params.clock_filter_iqd)
        .set{ refined_tree_ch }

    ancestral_sequences_ch = refined_tree_ch
        |combine(alignment_ch)
        |view()
        |map { arr -> tuple(params.inference_method, arr[0], arr[2]) }
        |ancestral_sequences
 
//    ancestral_sequences_ch |view 

//    ancestral_sequences(ancestral_sequences_ch)

    tree_ch_from_refined_tree = refined_tree_ch.map { it[0] }
    branch_ch_from_refined_tree = refined_tree_ch.map { it[1] }

    translate_ch = tree_ch_from_refined_tree
        |combine(ancestral_sequences.out)
        |combine(reference_ch)
        |translate

    traits_ch = tree_ch_from_refined_tree
        | combine(input_metadata_ch)
        | map { arr -> tuple(params.columns, arr[0], arr[1]) }
        | traits
    
    traits_ch.view()
    traits_from_traits_ch = traits_ch.map { it[0] }

    export_ch = tree_ch_from_refined_tree
        | combine(input_metadata_ch)
        | combine(branch_ch_from_refined_tree)
        | combine(traits_from_traits_ch)
        | combine(ancestral_sequences.out)
        | combine(translate.out)
        | combine(colors_ch)
        | combine(lat_longs_ch)
        | combine(auspice_config_ch)
        |export
}

