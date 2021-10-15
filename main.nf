#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams sample.bam [Options]

    Inputs Options:
    --input                     Input file: a csv file with paths to fastq files.
                                Example:
                                sample,fastq_1,fastq_2
                                SRA123456,path/to/SRA123456_1.fq.gz,path/to/SRA123456_2.fq.gz
                                SRA345678,path/to/SRA345678_1.fq.gz,path/to/SRA345678_2.fq.gz

    --genome                    Genome build version. Possible values: ${params.genomes.keySet().join(", ")}.

    Optional arguments:
    --adapters                  Fasta file contaning adapter sequences to be trimmed from the raw fastq files.

    --cleanup                   This option will enable nextflow work folder cleanup upon pipeline successfull completion.
                                All intermediate files from nexftlow processes' workdirs will be cleared, staging folder with staged
                                files will not be cleared.
                                If pipeline is completed with errors or interrupted cleanup will not be executed. Following successfull run
                                resumed from the failed run with --cleanup option enabled will only clear folders of processess created in
                                the latest run, it will not clear cached folders coming from previous pipleine runs.
                                (default: false)

    Resource Options:
    --max_cpus                  Maximum number of CPUs per task (int)
                                (default: $params.max_cpus)
    --max_memory                Maximum memory per task (memory unit)
                                (default: $params.max_memory)
    --max_time                  Maximum time per task (time unit)
                                (default: $params.max_time)

    See here for more info: https://github.com/zenomeplatform/nf-germline-snv
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}



/*
 * Functions
 */

// Function to ensure that resource requirements don't go beyond a maximum limit
def check_max(obj, type) {
  if (type == 'memory') {
    try {
      if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
        return params.max_memory as nextflow.util.MemoryUnit
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
      return obj
    }
  } else if (type == 'time') {
    try {
      if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
        return params.max_time as nextflow.util.Duration
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
      return obj
    }
  } else if (type == 'cpus') {
    try {
      return Math.min( obj, params.max_cpus as int )
    } catch (all) {
      println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
      return obj
    }
  }
}


// Simplified conditional known_sites parameters
known_sites         = params.known_sites ? params.known_sites : params.genomes[params.genome].containsKey("known_sites") ? params.genomes[params.genome].known_sites : null
known_sites_index   = params.known_sites_index ? params.known_sites_index : params.genomes[params.genome].containsKey("known_sites_index") ? params.genomes[params.genome].known_sites_index : null
known_sites_2       = params.known_sites_2 ? params.known_sites_2 : params.genomes[params.genome].containsKey("known_sites_2") ? params.genomes[params.genome].known_sites_2 : null
known_sites_2_index = params.known_sites_2_index ? params.known_sites_2_index : params.genomes[params.genome].containsKey("known_sites_2_index") ? params.genomes[params.genome].known_sites_2_index : null
known_sites_3       = params.known_sites_3 ? params.known_sites_3 : params.genomes[params.genome].containsKey("known_sites_3") ? params.genomes[params.genome].known_sites_3 : null
known_sites_3_index = params.known_sites_3_index ? params.known_sites_3_index : params.genomes[params.genome].containsKey("known_sites_3_index") ? params.genomes[params.genome].known_sites_3_index : null


// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Max Resources per task']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
summary['Input file']       = params.input
if (!params.fasta) summary['Genome build'] = params.genome
if (params.fasta) summary['Reference genome'] = params.fasta
if (params.fasta) summary['Reference genome index'] = params.fasta_fai
if (params.fasta) summary['Reference genome dict'] = params.fasta_dict
if (params.bwa) summary['BWA index'] = params.bwa
if (params.adapters) summary['Adapters'] = params.adapters
if (params.regions) summary['Target regions'] = params.regions
summary['Known sites']       = known_sites
summary['Known sites index'] = known_sites_index
if ( known_sites_2 ) {
  summary['Known sites 2']   = known_sites_2
  summary['Known sites 2 index'] = known_sites_2_index
}
if ( known_sites_3 ) {
  summary['Known sites 3']   = known_sites_3
  summary['Known sites 3 index'] = known_sites_3_index
}
if (params.bqsr_static_quantized_quals) summary['BQSR static quantized quals'] = params.bqsr_static_quantized_quals
if (params.cleanup) summary['Cleanup'] = "Cleanup is turned on"

log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"



// Define channels from repository files
projectDir = workflow.projectDir



/*
 * Check all important required inputs
 */

to_exit = false

if (!params.input) {to_exit=true; log.error "Please provide a csv file with paths to fastq files with --input [file] argument. See Readme for details."}

if (workflow.profile.contains('yandex')) {
  if (!params.accessKey) {to_exit=true; log.error "Please provide yandex accessKey to be used for yandex s3 bucket access with --accessKey '<key>' argument. Single quotes are important."}
  if (!params.secretKey) {to_exit=true; log.error "Please provide yandex secretKey to be used for yandex s3 bucket access with --secretKey '<key>' argument. Single quotes are important."}
}

// Check if genome parameter has a correct value
if (!params.genomes.keySet().contains(params.genome)) {
  to_exit=true; log.error "Reference data for genome \"${params.genome}\" is not available. Please choose one of: ${params.genomes.keySet().join(", ")}. \nAs an alternative you can provide your own reference files with parameters: \nnextflow run . [other parameters] --fasta file.fa --fasta_fai file.fa.fai --bwa file.fa.{amb,ann,bwt,pac,sa} \nNote: regex defition way for bwa has to be followed exactly as in example above, change only the basename to provide your own bwa index files."
}

// Check if at least one file with known variants is provided
if (!known_sites) {
  to_exit=true; log.error "Known genomic variants file has not been provided. It is a required input for GATK base recalibration step. \nPlease provide at least one known variants vcf file and its index with parameters --known_variants and --known_variants_index. \nAdditional variant files can be provided with --known_variants_2 and --known_variants_2_index; --known_variants_3 and --known_variants_3_index."
} else {
  if (!known_sites_index) {to_exit=true; log.error "Missing index file for genomic variants file \"${known_sites}\". It is a required input for GATK base recalibration step. \nPlease provide the index with parameter --known_variants_index."}
}

if (known_sites_2 && !known_sites_2_index) {to_exit=true; log.error "Missing index file for genomic variants file \"${known_sites_2}\". It is a required input for GATK base recalibration step. \nPlease provide the index with parameter --known_variants_2_index."}
if (known_sites_3 && !known_sites_3_index) {to_exit=true; log.error "Missing index file for genomic variants file \"${known_sites_3}\". It is a required input for GATK base recalibration step. \nPlease provide the index with parameter --known_variants_3_index."}


if (to_exit) exit 1, "One or more inputs are missing. Aborting."


/*
 * Define Channels from input
 */

Channel
    .fromPath(params.input)
    .ifEmpty { exit 1, "Cannot find input file : ${params.input}" }
    .splitCsv(skip:1)
    .map {sample_name, fastq_path_1, fastq_path_2, LB  ->

      if (params.metadata_from_file_name) {
          sample_name_parsed = file(fastq_path_1).simpleName.split(params.sample_name_format_delimeter)

          if (sample_name_parsed.length != params.sample_name_expected_item_number) {
            exit 1, "Error when parsing fastq file name. Each file name must have the following format: \n${params.sample_name_format}. \nExample         : ${params.sample_name_format_example} \nFailed file name: ${file(fastq_path_1).getName()} \n\nPlease follow the fastq file naming format, because pipeline extracts sample metadata (bio-type, seq-type, seq-machine, flowcell-ID, lane, and barcode) from the file name. \nYou can disable infering metadata from file name with parameter: \n--metadata_from_file_name false \nIn this case bio-type, seq-type and seq-machine will be left blank, and all samples will be processed as completely indepentdent samples."
          }

          if (params.double_check_sample_id) {
            sample_name_from_file_name = sample_name_parsed[0]
            if (sample_name != sample_name_from_file_name) {
              exit 1, "Error: sample name mismatch. \nSample name defined in input csv file sample_id column for file ${file(fastq_path_1).getName()} does not match the one specified in the file name. \nFrom sample_id column: ${sample_name} \nFrom file name       : ${sample_name_from_file_name} \n\nPlease provide correct sample and file name, or disable the sample name checking with parameter: \n--double_check_sample_id false"
            }
          }

          bio_type    = sample_name_parsed[1]
          seq_type    = sample_name_parsed[2]
          seq_machine = sample_name_parsed[3]
          flowcell_id = sample_name_parsed[4]
          lane        = sample_name_parsed[5]
          barcode     = sample_name_parsed[6]
      }

      if (!params.metadata_from_file_name) {
          bio_type    = "not_provided"
          seq_type    = "not_provided"
          seq_machine = "IL"
          flowcell_id = sample_name
          lane        = sample_name
          barcode     = sample_name
      }

      read_group_LB = LB

      [ sample_name, file(fastq_path_1), file(fastq_path_2), bio_type, seq_type, seq_machine, flowcell_id, lane, barcode, read_group_LB ]
    }
    .set { ch_input_fastq }

ch_input_fastq.into {
  ch_input_fastq_for_qc;
  ch_input_fastq_to_trim
}

refgenome = params.fasta ? params.fasta : params.genomes[params.genome].fasta
refgenome_index = params.fasta ? params.fasta_fai : params.genomes[params.genome].fasta_fai
refgenome_dict = params.fasta ? params.fasta_dict : params.genomes[params.genome].fasta_dict

ch_refgenome = Channel.value(file(refgenome))
ch_refgenome_index = Channel.value(file(refgenome_index))
ch_refgenome_dict = Channel.value(file(refgenome_dict))

bwa = params.bwa ? params.bwa : params.genomes[params.genome].bwa
ch_bwa = Channel.fromFilePairs(bwa, size: 5, flat: true)

ch_known_sites_1 = Channel.value(file(known_sites))
ch_known_sites_1_index = Channel.value(file(known_sites_index))

// Optional inputs
ch_adapters = params.adapters ? Channel.value(file(params.adapters)) : "null"

if (known_sites_2) {
  ch_known_sites_2 = Channel.value(file(known_sites_2))
  ch_known_sites_2_index = Channel.value(file(known_sites_2_index))
} else {
  ch_known_sites_2 = "null"
  ch_known_sites_2_index = "null"
}

if (known_sites_3) {
  ch_known_sites_3 = Channel.value(file(known_sites_3))
  ch_known_sites_3_index = Channel.value(file(known_sites_3_index))
} else {
  ch_known_sites_3 = "null"
  ch_known_sites_3_index = "null"
}

ch_regions = params.regions ? Channel.value(file(params.regions)) : "null"


/*
 * Processes
 */

process fastqc_raw {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/fastqc/raw/", mode: 'copy'

    input:
    set val(sample_name),
        file(fastq_1),
        file(fastq_2),
        val(bio_type),
        val(seq_type),
        val(seq_machine),
        val(flowcell_id),
        val(lane),
        val(barcode),
        val(read_group_LB) from ch_input_fastq_for_qc

    output:
    set val(sample_name), file("fastqc_${sample_name}_raw_logs") into ch_fastq_qc_raw


    script:
    """
    mkdir fastqc_${sample_name}_raw_logs

    fastqc --outdir fastqc_${sample_name}_raw_logs \
        --format fastq \
        --quiet \
        --threads ${task.cpus} \
        ${fastq_1} ${fastq_2}
    """
  }


process trim_fastqc {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/trimmed_reads/", mode: 'copy'

    input:
    set val(sample_name),
        file(fastq_1),
        file(fastq_2),
        val(bio_type),
        val(seq_type),
        val(seq_machine),
        val(flowcell_id),
        val(lane),
        val(barcode),
        val(read_group_LB) from ch_input_fastq_to_trim
    each file(adapters) from ch_adapters

    output:
    set val(sample_name),
        file("${fastq_1.simpleName}_trim.fastq.gz"),
        file("${fastq_2.simpleName}_trim.fastq.gz"),
        val(bio_type),
        val(seq_type),
        val(seq_machine),
        val(flowcell_id),
        val(lane),
        val(barcode),
        val(read_group_LB) into (ch_fastq_trimmed_to_map, ch_fastq_trimmed_for_qc)
    set val(sample_name), file("flexbar_${sample_name}.log") into ch_trimming_report

    script:
    // adapters were not provide so for now are optional
    adapters_param = params.adapters ? "--adapters ${adapters}" : ""
    """
    flexbar --adapter-min-overlap ${params.adapter_min_overlap} \
        --adapter-trim-end ${params.adapter_trim_end} \
        --pre-trim-left ${params.pre_trim_left} \
        --max-uncalled ${params.max_uncalled} \
        --min-read-length ${params.min_read_length} \
        --threads ${task.cpus} \
        --zip-output GZ \
        --reads ${fastq_1} \
        --reads2 ${fastq_2} \
        --output-reads ${fastq_1.simpleName}_trim.fastq \
        --output-reads2 ${fastq_2.simpleName}_trim.fastq \
        --output-log flexbar_${sample_name}.log \
        $adapters_param

        # no .gz in the end of output files is important, its added automatically by flexbar because of --zip-output GZ option.
    """
  }


process fastqc_trimmed {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/fastqc/trimmed/", mode: 'copy'

    input:
    set val(sample_name),
        file(fastq_1),
        file(fastq_2),
        val(bio_type),
        val(seq_type),
        val(seq_machine),
        val(flowcell_id),
        val(lane),
        val(barcode),
        val(read_group_LB) from ch_fastq_trimmed_for_qc

    output:
    set val(sample_name), file("fastqc_${sample_name}_trimmed_logs") into ch_fastq_qc_trimmed

    script:
    """
    mkdir fastqc_${sample_name}_trimmed_logs

    fastqc --outdir fastqc_${sample_name}_trimmed_logs \
        --format fastq \
        --quiet \
        --threads ${task.cpus} \
        ${fastq_1} ${fastq_2}
    """
  }


ch_fastq_qc_raw
    .join(ch_trimming_report, by: 0)
    .join(ch_fastq_qc_trimmed, by: 0)
    .into { ch_prealignment_multiqc_files_by_sample; ch_prealignment_multiqc_files_all }


if (params.multiqc_prealignment_by_sample) {
  process multiqc_prealignment_report_by_sample {
      tag "$sample_name"
      label 'low_memory'
      publishDir "${params.outdir}/multiqc_prealignment_report/${sample_name}/", mode: 'copy'

      input:
      set val(sample_name), file(fastqc_raw_dir), file(trimming_log), file(fastqc_trimmed_dir) from ch_prealignment_multiqc_files_by_sample

      output:
      file("multiqc_report.html")

      script:
      """
      multiqc .
      """
    }
}


if (params.multiqc_prealignment_all) {
  process multiqc_prealignment_report_all {
      label 'low_memory'
      publishDir "${params.outdir}/multiqc_prealignment_report/", mode: 'copy'

      input:
      file("*") from ch_prealignment_multiqc_files_all.collect()

      output:
      file("multiqc_report.html")

      script:
      """
      multiqc .
      """
    }
}


process map_reads {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${sample_name}/align/", mode: 'copy'

    input:
    set val(sample_name),
        file(fastq_1),
        file(fastq_2),
        val(bio_type),
        val(seq_type),
        val(seq_machine),
        val(flowcell_id),
        val(lane),
        val(barcode),
        val(read_group_LB) from ch_fastq_trimmed_to_map
    each file(bwa_indexes) from ch_bwa

    output:
    set val(sample_name), file("${fastq_1.simpleName}_L${read_group_LB}.bam") into ch_mapped_reads

    script:
    seq_type_modified = sample_name[2]

    if (params.seq_machine_catalog.keySet().contains(seq_machine)) {
      seq_machine_modified = params.seq_machine_catalog[seq_machine].full_name
    } else {
      log.warn "Sample ${fastq_1.name} has unknown sequening machine type \"${seq_machine}\". \nKnow machines in catalog: ${params.seq_machine_catalog.keySet().join(", ")}. \nSample will receive \"Unknown_seq_machine\" value."
      seq_machine_modified = "Unknown_seq_machine"
    }

    """
    bwa mem -t ${task.cpus} \
        -Y \
        -R "@RG\\tID:${seq_type_modified}\\tPL:${seq_machine_modified}\\tPU:v${flowcell_id}.${lane}.${barcode}\\tLB:exome_lib${read_group_LB}\\tSM:${sample_name}" \
        ${bwa_indexes[1].baseName} \
        ${fastq_1} \
        ${fastq_2} \
      | samtools view -bS -@${task.cpus} - > ${fastq_1.simpleName}_L${read_group_LB}.bam;
    """
}


ch_mapped_reads_grouped_by_sample = ch_mapped_reads.groupTuple(by: 0).view()


process merge_bams_by_sample {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${sample_name}/align/", mode: 'copy'

    input:
    set val(sample_name), file(bam_files) from ch_mapped_reads_grouped_by_sample

    output:
    set val(sample_name), file("${sample_name}.merged.bam") into ch_mapped_reads_merged_by_sample

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_name}.merged.bam ${bam_files}
    """
}


process sort_bams_by_name {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${sample_name}/align/", mode: 'copy'

    input:
    set val(sample_name), file(merged_bam) from ch_mapped_reads_merged_by_sample

    output:
    set val(sample_name), file("${sample_name}.namesorted.bam") into ch_mapped_reads_sorted_by_name

    script:
    """
    picard SortSam \
        I=${merged_bam} \
        O=${sample_name}.namesorted.bam \
        SO=queryname
    """
}


process mark_duplicates  {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${sample_name}/align/", mode: 'copy'

    input:
    set val(sample_name), file(sorted_bam) from ch_mapped_reads_sorted_by_name

    output:
    set val(sample_name), file("${sample_name}.namesorted_mrkdup.bam") into ch_mapped_reads_mrkdup
    file("${sample_name}_mrkdup_metrics.txt") into ch_ch_mapped_reads_mrkdup_metrics

    script:
    """
    picard MarkDuplicates \
        I=${sorted_bam} \
        O=${sample_name}.namesorted_mrkdup.bam \
        ASSUME_SORT_ORDER=queryname \
        METRICS_FILE=${sample_name}_mrkdup_metrics.txt \
        QUIET=true \
        COMPRESSION_LEVEL=0 \
        VALIDATION_STRINGENCY=LENIENT
    """
}


process sort_bams_by_coord  {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${sample_name}/align/", mode: 'copy'

    input:
    set val(sample_name), file(bam_mrkdupl) from ch_mapped_reads_mrkdup

    output:
    set val(sample_name), file("${sample_name}.sorted_mrkdup.bam"), file("${sample_name}.sorted_mrkdup.bai") into ch_mapped_reads_mrkdup_sorted

    script:
    """
    picard SortSam \
        I=${bam_mrkdupl} \
        O=${sample_name}.sorted_mrkdup.bam \
        SO=coordinate

    picard BuildBamIndex \
        I=${sample_name}.sorted_mrkdup.bam
    """
}

process calculate_BQSR  {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${sample_name}/align/", pattern: "*.table", mode: 'copy'

    input:
    set val(sample_name), file(bam), file(bam_index) from ch_mapped_reads_mrkdup_sorted
    file(fasta) from ch_refgenome
    file(fasta_fai) from ch_refgenome_index
    file(fasta_dict) from ch_refgenome_dict
    file(known_sites_1) from ch_known_sites_1
    file(known_sites_1_index) from ch_known_sites_1_index
    each file(known_sites_2) from ch_known_sites_2
    each file(known_sites_2_index) from ch_known_sites_2_index
    each file(known_sites_3) from ch_known_sites_3
    each file(known_sites_3_index) from ch_known_sites_3_index
    each file(regions) from ch_regions

    output:
    set val(sample_name), file(bam), file(bam_index), file("${sample_name}.sorted_mrkdup_bqsr.table") into ch_mapped_reads_with_BQSR

    script:
    optional_known_sites_2_arg = params.known_sites_2 ? "--known-sites ${known_sites_2}" : ""
    optional_known_sites_3_arg = params.known_sites_3 ? "--known-sites ${known_sites_3}" : ""
    optional_regions_file      = params.regions ? "-L ${regions}" : ""
    """
    gatk BaseRecalibrator \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample_name}.sorted_mrkdup_bqsr.table \
        --known-sites ${known_sites_1} \
        ${optional_known_sites_2_arg} \
        ${optional_known_sites_3_arg} \
        --preserve-qscores-less-than ${params.bqsr_preserve_qscores_less_than} \
        --disable-bam-index-caching \
        ${optional_regions_file}
    """
}

process apply_BQSR  {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${sample_name}/align/", pattern: "*.table", mode: 'copy'

    input:
    set val(sample_name), file(bam), file(bam_index), file(bqsr_table) from ch_mapped_reads_with_BQSR
    file(fasta) from ch_refgenome
    file(fasta_fai) from ch_refgenome_index
    file(fasta_dict) from ch_refgenome_dict

    output:
    set val(sample_name), file("${sample_name}.sorted_mrkdup_bqsr.bam"), file("${sample_name}.sorted_mrkdup_bqsr.bai") into ch_recalibrated_mapped_reads

    script:
    arg_static_quantized_quals = ""
    if (params.bqsr_static_quantized_quals) {
      static_quantized_quals_array = params.bqsr_static_quantized_quals.toString().split(',')
      for(qual in static_quantized_quals_array) {
         arg_static_quantized_quals += "--static-quantized-quals ${qual} "
      }
    }
    """
    gatk ApplyBQSR \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample_name}.sorted_mrkdup_bqsr.bam \
        --bqsr-recal-file ${bqsr_table} \
        --preserve-qscores-less-than ${params.bqsr_preserve_qscores_less_than} \
        ${arg_static_quantized_quals}

    picard BuildBamIndex \
        I=${sample_name}.sorted_mrkdup_bqsr.bam
    """
}


/*
 * Workflow completion
 */

workflow.onComplete {

    c_green = "\033[0;32m";
    c_purple = "\033[0;35m";
    c_red = "\033[0;31m";
    c_reset = "\033[0m";

    if (workflow.success) {
        log.info "-${c_purple}[splicing-pipelines-nf]${c_green} Pipeline completed successfully${c_reset}-"
        if (params.cleanup) {
          log.info "-${c_purple}[splicing-pipelines-nf]${c_green} Cleanup: Working directory cleared from intermediate files generated with current run: '${workflow.workDir}'  ${c_reset}-"
        }
    } else { // To be shown requires errorStrategy = 'finish'
        log.info "-${c_purple}[splicing-pipelines-nf]${c_red} Pipeline completed with errors${c_reset}-"
        if (params.cleanup) {
          log.info "-${c_purple}[splicing-pipelines-nf]${c_red} Cleanup: Working directory was not cleared from intermediate files due to pipeline errors. You can re-use them with -resume option.${c_reset}-"
        }
    }
}
