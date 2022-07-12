
include { GATK4_CREATESEQUENCEDICTIONARY }              from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER }                       from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
// include { GET_META }                                    from '../../modules/local/getmeta/main'
// include { GET_PATH as GET_PATH_ARRIBA_FAIL }            from '../../modules/local/getpath/main'
// include { SAMTOOLS_SORT as SAMTOOLS_SORT_FOR_ARRIBA }   from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_GATK}    from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_FAIDX }                              from '../../modules/nf-core/modules/samtools/faidx/main'


workflow GATK_WORKFLOW {
    take:
        reads
        ch_fasta
        bam_indexed

    main:
        ch_versions = Channel.empty()
        ch_dict = Channel.empty()

        if ((params.gatk || params.all) && !params.fusioninspector_only) {

            GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
            ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
            ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

            SAMTOOLS_FAIDX(ch_fasta)
            ch_fai = SAMTOOLS_FAIDX.out.fai
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

            SAMTOOLS_INDEX(ch_fasta)
            ch_fai = SAMTOOLS_INDEX.out.fai
            ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)


            GATK4_HAPLOTYPECALLER()

            STAR_FOR_ARRIBA( reads, ch_starindex_ref, ch_gtf, params.star_ignore_sjdbgtf, '', params.seq_center ?: '')
            ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions)

        }
        else {
            ch_arriba_fusions       = GET_META(reads, ch_dummy_file)
            ch_arriba_fusion_fail   = ch_dummy_file
            ch_arriba_visualisation = ch_dummy_file
        }

    emit:
        fusions         = ch_arriba_fusions
        fusions_fail    = ch_arriba_fusion_fail
        versions        = ch_versions.ifEmpty(null)
        pdf             = ch_arriba_visualisation
    }

