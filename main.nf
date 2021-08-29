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
if (params.adapters) summary['Adapters'] = params.adapters
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

if (to_exit) exit 1, "One or more inputs are missing. Termingating pipeline."


/*
 * Define Channels from input
 */

Channel
    .fromPath(params.input)
    .ifEmpty { exit 1, "Cannot find input file : ${params.input}" }
    .splitCsv(skip:1)
    .map {sample_name, fastq_path_1, fastq_path_2  -> [ sample_name, file(fastq_path_1), file(fastq_path_2) ] }
    .set { ch_input_fastq }

ch_input_fastq.into {
  ch_input_fastq_for_qc;
  ch_input_fastq_to_trim
}

// Optional inputs
ch_adapters = params.adapters ? Channel.value(file(params.adapters)) : "null"



/*
 * Processes
 */

process fastqc_raw {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/fastqc/raw/", mode: 'copy'

    input:
    set val(sample_name), file(fastq_1), file(fastq_2) from ch_input_fastq_for_qc

    output:
    file("fastqc_${sample_name}_raw_logs") into ch_fastq_qc_raw


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
    set val(sample_name), file(fastq_1), file(fastq_2) from ch_input_fastq_to_trim
    each file(adapters) from ch_adapters

    output:
    set val(sample_name), file("${fastq_1.simpleName}_trimmed.fastq.gz"), file("${fastq_2.simpleName}_trimmed.fastq.gz") into (ch_fastq_trimmed_to_map, ch_fastq_trimmed_for_qc)
    file("flexbar_${sample_name}.log") into ch_trimming_report

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
        --output-reads ${fastq_1.simpleName}_trimmed.fastq \
        --output-reads2 ${fastq_2.simpleName}_trimmed.fastq \
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
    set val(sample_name), file(fastq_1), file(fastq_2) from ch_fastq_trimmed_for_qc

    output:
    file("fastqc_${sample_name}_trimmed_logs") into ch_fastq_qc_trimmed

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


process multiqc_prealignment_report {
    label 'low_memory'
    publishDir "${params.outdir}/multiqc_prealignment_report/", mode: 'copy'

    input:
    file(fastqc_raw_dir) from ch_fastq_qc_raw.collect()
    file(trimming_log) from ch_trimming_report.collect()
    file(fastqc_trimmed_dir) from ch_fastq_qc_trimmed.collect()

    output:
    file("multiqc_report.html")

    script:
    """
    multiqc .
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
