include {FASTQC} from './modules/fastqc'
include {STAR_INDEX} from './modules/star_index'
include {STAR_ALIGN} from './modules/star_align'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {PICARD_CREATE_DICT} from './modules/picard_create_dict'
include {GTF_TO_REFFLAT} from './modules/gtf_to_refflat'
include {GTF_TO_RRNA_BED} from './modules/gtf_to_rrna_bed'
include {BED_TO_INTERVAL_LIST} from './modules/bed_to_interval_list'
include {PICARD_COLLECT_RNASEQ_METRICS} from './modules/picard_collect_rnaseq_metrics'
include {MULTIQC} from './modules/multiqc'
include {SAMTOOLS_SORT} from './modules/samtools_sort'
include {BAM_INDEX} from './modules/bam_index'
include {MERGE_BAMS} from './modules/merge_bams'
include {SAMTOOLS_PILEUP} from './modules/samtools_pileup'
include {VARSCAN_TUMOR} from './modules/varscan_somatic_tumor'
include {VARSCAN_NONTUMOR} from './modules/varscan_somatic_nontumor'
include {FILTER_VARIANTS} from './modules/filter_variants'
include {BUILD_SNPEFF_DB} from './modules/build_snpeff_db'
include {ANNOTATE_VARIANTS} from './modules/annotate_variants'
include {PROCESS_VCFS} from './modules/process_vcfs'
include {CREATE_OUTPUT} from './modules/create_output'
include {CSV_TO_XLSX} from './modules/csv_to_xlsx'

workflow {
    Channel.fromPath(params.samplesheet)
    | splitCsv( header: true )
    | map{ row -> tuple(row.name, file(row.read1), file(row.read2)) }
    | set { reads_ch }

    FASTQC(reads_ch)

    if ( params.star_index && file(params.star_index).exists() ) {
        star_index_ch = Channel.value( file(params.star_index) )
    } else {
        STAR_INDEX(
            Channel.value( file(params.gtf) ),
            Channel.value( file(params.ref_genome) )
        )
        star_index_ch = STAR_INDEX.out.index
    }

    STAR_ALIGN(reads_ch, star_index_ch)

    SAMTOOLS_FLAGSTAT(STAR_ALIGN.out.bam)
    
    if ( params.fa_dict && file(params.fa_dict).exists() ) {
        fa_dict_ch = Channel.value( file(params.fa_dict) )
    } else {
        PICARD_CREATE_DICT( Channel.value( file(params.ref_genome) ) )
        fa_dict_ch = PICARD_CREATE_DICT.out.dict
    }

    GTF_TO_REFFLAT(params.gtf)

    // This is currently blank because there is no rRNA inside the mm9 GTF file, ask max about this once I have a better idea
    GTF_TO_RRNA_BED(params.gtf)
    
    BED_TO_INTERVAL_LIST(GTF_TO_RRNA_BED.out.bed, fa_dict_ch)

    // Picard RNA metrics per sample
    PICARD_COLLECT_RNASEQ_METRICS(
      STAR_ALIGN.out.bam,
      GTF_TO_REFFLAT.out.refflat,
      BED_TO_INTERVAL_LIST.out.intervals,
      params.ref_genome
    )
    
    multiqc_ch = FASTQC.out.zip
        .map { it[1] }
        .mix(SAMTOOLS_FLAGSTAT.out)
        .mix(PICARD_COLLECT_RNASEQ_METRICS.out.metrics)
        .collect()

    //multiqc_ch.view()
    
    MULTIQC(multiqc_ch)

    SAMTOOLS_SORT(STAR_ALIGN.out.bam)

    BAM_INDEX(SAMTOOLS_SORT.out)

    BAM_INDEX.out
        .branch {
            control: it[0].startsWith('Control')
            exp: !it[0].startsWith('Control')
        }
        .set { branched_channels }

    control_ch = branched_channels.control
    exp_ch = branched_channels.exp

    control_ch
        .map { sample, bam, bai -> ['merged_controls', bam, bai] }
        .groupTuple()
        .set { all_controls_ch }

    MERGE_BAMS(all_controls_ch)

    exp_ch
        .combine(MERGE_BAMS.out)
        .set { paired_ch }
    
    SAMTOOLS_PILEUP(paired_ch, params.ref_genome)

    SAMTOOLS_PILEUP.out
        .branch {
            tumor: it[2].startsWith('Tumor')
            exp_nontumor: !it[2].startsWith('Tumor')
        }
        .set { paired_branched_channels }

    exp_nontumor_ch = paired_branched_channels.exp_nontumor
    tumor_ch = paired_branched_channels.tumor

    VARSCAN_NONTUMOR(exp_nontumor_ch, params.control_cnt)

    VARSCAN_TUMOR(tumor_ch, params.control_cnt)

    varscan_combined_ch = VARSCAN_NONTUMOR.out.snp.mix(VARSCAN_TUMOR.out.snp)
        .mix(VARSCAN_NONTUMOR.out.indel)
        .mix(VARSCAN_TUMOR.out.indel)

    FILTER_VARIANTS(varscan_combined_ch)

    BUILD_SNPEFF_DB(params.ref_genome, params.gtf)

    BUILD_SNPEFF_DB.out.view()

    ANNOTATE_VARIANTS(FILTER_VARIANTS.out, BUILD_SNPEFF_DB.out)
    
    def lncrna_file = (params.lncRNAs_ref && file(params.lncRNAs_ref).exists())
        ? file(params.lncRNAs_ref)
        : file("$projectDir/bin/NO_FILE_LNCRNA")

    def known_snps_file = (params.known_snps_vcf && file(params.known_snps_vcf).exists())
        ? file(params.known_snps_vcf)
        : file("$projectDir/bin/NO_FILE_KNOWNSNPS")

    lncrna_ch    = Channel.value(lncrna_file)
    known_snps_ch = Channel.value(known_snps_file)

    PROCESS_VCFS(ANNOTATE_VARIANTS.out, lncrna_ch, known_snps_ch)

    def groupParts = { String id ->
        // id examples: control_vs_aHepGHRkd_Stat5bCA4_indel  OR  control_vs_aHepGHRkd3_snp
        def core = id.replaceFirst(/_(snp|indel)$/, '')   // drop _snp/_indel
        def m = (core =~ /^(.*?)(\d+)$/)                 // trailing digits at end
        if( m.matches() ) {
            return [ m[0][1], (m[0][2] as Integer) ]     // [groupPrefix, replicateNum]
        } else {
            return [ core, Integer.MAX_VALUE ]           // no replicate number
        }
    }

    def grouped = PROCESS_VCFS.out
        .map { id, csv ->
            def (grp, num) = groupParts(id)
            tuple(grp, num, id, csv)   // (group, num, full_id, file)
        }

    snp_by_group = grouped
        .filter { grp, num, id, csv -> id.endsWith('snp') }
        .map { grp, num, id, csv -> tuple(grp, num, csv) }      // (group, num, file)
        .groupTuple()
        .map { grp, nums, files ->
            def ordered = [nums, files].transpose()
                .sort { it[0] }
                .collect { it[1] }                              // List<Path>
            tuple(grp, ordered)
        }

    indel_by_group = grouped
        .filter { grp, num, id, csv -> id.endsWith('indel') }
        .map { grp, num, id, csv -> tuple(grp, num, csv) }
        .groupTuple()
        .map { grp, nums, files ->
            def ordered = [nums, files].transpose()
                .sort { it[0] }
                .collect { it[1] }
            tuple(grp, ordered)
        }

    paired = indel_by_group
        .join(snp_by_group)                         
        .map { grp, indels, snps -> tuple(grp, indels, snps) }

    CREATE_OUTPUT(paired)

    individuals_excel_ch = paired.map { id, file1, file2 -> [file1, file2] }.flatten()

    CSV_TO_XLSX(individuals_excel_ch)
}
