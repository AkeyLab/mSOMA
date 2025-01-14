nextflow.enable.dsl=2

// Required inputs
params.sampleSheet = file_exists(param_set(params.sampleSheet,"--sampleSheet"))
params.outputDir = file(param_set(params.outputDir, "--outputDir"))

// Default parameters
// can overwride when running nextflow with --min-MQ 30 for example
params.minMQ = 20
params.requireFlags = 2
params.excludeFlags = 3844
params.maxIndel = 0
params.mismatchFrac = 0.1
params.softclipFrac = 0.5
params.seqLength = 75
params.ntrim = 5
params.minBQ = 20
params.maxAltAllele = 1
params.minDepth = 10
params.help = false

// Provide a help message if someone runs `nextflow run --help script.nf`
if (params.help) {
    help = """mSOMA somatic mutation caller
             |
             |Required arguments:
             |  --sample-sheet    Path to a tab-delimited file with the following columns:
             |                    Sample, Bam, Bai, Bed, Fasta
             |
             |  --output-dir          Path to the output directory
             |
             |Optional arguments:
             |  Refer to the mSOMA documentation for more information.
           """.stripMargin()

    println(help)
    exit(0)
}

// Write info to stdout about the current run
log.info """\
         mSOMA
         ==========================
         input from   : ${params.sampleSheet}
         output to    : ${params.outputDir}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()

// Helper functions to validate inputs and provide helpful error messages
def param_set(param, name) {
    if (!param) {
        System.err.println("ERROR: required input "+name+" is not set")
        System.exit(1)
    }
    return param
}
def file_exists(f) {
    f = file(f)
    if (!f.exists()) {
        System.err.println("ERROR in validate_path: can't find file: "+f)
        System.exit(1)
    }
    return f
}

process countReads {
    publishDir params.outputDir, mode: 'copy'
    cache 'lenient'

    input:
        tuple val(s), file(bam), file(bai), file(bed), file(fasta)
    output:
        tuple val(s), file("${s}.counts.gz")
    """
    msoma count \
        --bed $bed \
        --fasta $fasta \
        --min-MQ $params.minMQ \
        --require-flags $params.requireFlags \
        --exclude-flags $params.excludeFlags \
        --max-indel $params.maxIndel \
        --mismatch-frac $params.mismatchFrac \
        --softclip-frac $params.softclipFrac \
        --seq-length $params.seqLength \
        --ntrim $params.ntrim \
        --min-BQ $params.minBQ \
        --min-depth $params.minDepth \
        --max-alt-allele $params.maxAltAllele \
        --output "${s}.counts.gz" \
        $bam
    """
}

process mle {
    label 'high_mem'

    publishDir params.outputDir, mode: 'copy'
    cache 'lenient'

    input:
        tuple val(s), file(counts)
    output:
        tuple file("${s}.pvals"), file("${s}.ab")
    """
    msoma mle \
        --output "${s}.pvals" \
        --ab "${s}.ab" \
        --min-depth $params.minDepth \
        $counts
    """
}

workflow {
    //Read in the samplesheet table
    rows = Channel.fromPath(params.sampleSheet) \
    | splitCsv(header:true, sep:"\t") \

    //Generate count files
    counts = rows \
    | map { r -> tuple(
         r.Sample,
         file(r.Bam),
         file(r.Bai),
         file(r.Bed),
         file(r.Fasta),
         ) } \
    | countReads

    //Perform mle a/b estimation and pval calculations
    counts | mle
}
