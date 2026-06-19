/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                              } from '../modules/nf-core/fastqc/main'
include { FASTP                               } from '../modules/nf-core/fastp/main'
include { BOWTIE2_BUILD                       } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                       } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_FAIDX                      } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                      } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                       } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW                       } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_RAW  } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_TRIM } from '../modules/nf-core/samtools/stats/main'
include { MOSDEPTH as MOSDEPTH_RAW            } from '../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_TRIM           } from '../modules/nf-core/mosdepth/main'
include { IVAR_TRIM                           } from '../modules/nf-core/ivar/trim/main'
include { IVAR_VARIANTS                       } from '../modules/nf-core/ivar/variants/main'
include { IVAR_CONSENSUS                      } from '../modules/nf-core/ivar/consensus/main'
include { LOFREQ_INDELQUAL                    } from '../modules/nf-core/lofreq/indelqual/main'
include { LOFREQ_CALL                         } from '../modules/nf-core/lofreq/call/main'
include { LOFREQ_FILTER                       } from '../modules/nf-core/lofreq/filter/main'
include { FIND_CONCATENATE                    } from '../modules/nf-core/find/concatenate/main'
include { MAFFT_ALIGN                         } from '../modules/nf-core/mafft/align/main'
include { IQTREE                              } from '../modules/nf-core/iqtree/main'
include { MULTIQC                             } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_arbor_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ARBOR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SHARED REFERENCE — one reference for the whole run, as value channels
    //
    def ch_ref       = channel.value([ [id:'reference'], file(params.reference, checkIfExists:true) ])
    def ch_ref_fasta = file(params.reference, checkIfExists:true)              // bare path (ivar/lofreq)
    def ch_primer_bed = file(params.primer_bed, checkIfExists:true)           // bare path (ivar/trim)
    def ch_context   = params.context_fasta
        ? channel.value([ [id:'context'], file(params.context_fasta, checkIfExists:true) ])
        : channel.value([ [], [] ])
    def ch_adapter   = params.adapter_fasta ? file(params.adapter_fasta, checkIfExists:true) : []

    BOWTIE2_BUILD(ch_ref)
    def ch_index = BOWTIE2_BUILD.out.index.first()

    // REVIEW: faidx fai feeds samtools stats/view/sort (tuple) and ivar variants (bare path)
    SAMTOOLS_FAIDX(ch_ref.map { m, f -> [ m, f, [] ] }, false)
    def ch_fai_tuple = SAMTOOLS_FAIDX.out.fai.map { m, fai -> [ m, ch_ref_fasta, fai ] }.first()
    def ch_fai_bare  = SAMTOOLS_FAIDX.out.fai.map { _m, fai -> fai }.first()

    //
    // READS QC + TRIM
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_samplesheet)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _m, f -> f })
    }
    // fastp takes [meta, reads, adapter_fasta] (adapter [] = auto-detect Illumina adapters)
    FASTP(ch_samplesheet.map { meta, reads -> [ meta, reads, ch_adapter ] }, false, false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ _m, f -> f })

    //
    // MAP (local) -> sorted BAM -> index for .bai
    // REVIEW: bowtie2 sort_bam=true emits .csi; ivar/trim needs .bai -> SAMTOOLS_INDEX
    //
    BOWTIE2_ALIGN(FASTP.out.reads, ch_index, ch_ref, false, true)
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.map{ _m, f -> f })
    SAMTOOLS_INDEX(BOWTIE2_ALIGN.out.bam)
    def ch_mapped = BOWTIE2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.index)      // [meta, bam, bai]

    // pre-trim mapping QC
    SAMTOOLS_STATS_RAW(ch_mapped, ch_fai_tuple)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_RAW.out.stats.map{ _m, f -> f })
    MOSDEPTH_RAW(ch_mapped.map { m, b, i -> [ m, b, i, [] ] }, ch_ref, [])
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_RAW.out.summary_txt.map{ _m, f -> f })

    //
    // TILING PRIMER TRIM -> re-sort + index
    //
    IVAR_TRIM(ch_mapped, ch_primer_bed)
    SAMTOOLS_SORT(IVAR_TRIM.out.bam, [[],[],[]], 'bai')
    def ch_final = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.index)        // [meta, bam, bai]

    // post-trim QC
    SAMTOOLS_STATS_TRIM(ch_final, ch_fai_tuple)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_TRIM.out.stats.map{ _m, f -> f })
    MOSDEPTH_TRIM(ch_final.map { m, b, i -> [ m, b, i, [] ] }, ch_ref, [])
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_TRIM.out.summary_txt.map{ _m, f -> f })

    //
    // VARIANTS (whole multi-segment reference)
    //
    if (!params.skip_ivar_variants) {
        // ivar/variants: [meta, bam], path fasta, path fai, path gff, val save_mpileup
        IVAR_VARIANTS(ch_final.map { m, b, _i -> [ m, b ] }, ch_ref_fasta, ch_fai_bare, [], false)
        ch_multiqc_files = ch_multiqc_files.mix(IVAR_VARIANTS.out.tsv.map{ _m, f -> f })
    }
    if (!params.skip_lofreq) {
        // REVIEW: lofreq indelqual pre-processing before call (recommended for Illumina)
        LOFREQ_INDELQUAL(ch_final.map { m, b, _i -> [ m, b ] }, ch_ref)
        LOFREQ_CALL(LOFREQ_INDELQUAL.out.bam.map { m, b -> [ m, b, [] ] }, ch_ref_fasta)
        LOFREQ_FILTER(LOFREQ_CALL.out.vcf)
    }

    //
    // PER-SEGMENT CONSENSUS + PHYLOGENY (the fan-out/fan-in — riskiest wiring)
    //
    if (!params.skip_phylogeny) {
        def ch_seg = channel.fromList(params.segments.tokenize(','))         // e.g. S, M, L
        // REVIEW: segment names MUST match reference FASTA record IDs exactly
        def ch_seg_bam = ch_final.combine(ch_seg)
            .map { meta, bam, bai, seg ->
                [ meta + [ id:"${meta.id}_${seg}", sample:meta.id, segment:seg ], bam, bai ]
            }
        // split BAM to one segment (region passed via ext.args2 = meta.segment)
        SAMTOOLS_VIEW(ch_seg_bam, ch_fai_tuple, [[],[]], [[],[]], '')
        // per-segment consensus (iVar is single-reference; the BAM is now one segment)
        IVAR_CONSENSUS(SAMTOOLS_VIEW.out.bam, ch_ref_fasta, false)
        // REVIEW: regroup consensus across samples BY segment -> one FASTA per segment
        def ch_seg_consensus = IVAR_CONSENSUS.out.fasta
            .map { meta, fa -> [ meta.segment, fa ] }
            .groupTuple()
            .map { seg, fas -> [ [ id:seg ], fas ] }
        FIND_CONCATENATE(ch_seg_consensus)
        // MSA per segment, folding in optional external reference strains via MAFFT --add
        MAFFT_ALIGN(FIND_CONCATENATE.out.file_out, ch_context, [[],[]], [[],[]], [[],[]], [[],[]], false)
        // ML tree per segment
        IQTREE(MAFFT_ALIGN.out.fas.map { meta, aln -> [ meta, aln, [] ] },
            [], [], [], [], [], [], [], [], [], [], [], [])
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'arbor_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'arbor'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
