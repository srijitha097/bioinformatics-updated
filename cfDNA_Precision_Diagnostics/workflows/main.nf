nextflow.enable.dsl = 2

process INDEX_REFERENCE {

    tag 'GRCh38 targeted reference'

    input:
    path reference

    output:
    path 'indexed_reference'

    script:
    """
    mkdir indexed_reference

    cp ${reference} indexed_reference/reference.fa

    bwa index indexed_reference/reference.fa

    samtools faidx indexed_reference/reference.fa
    """
}


process FASTQC {

    tag "$sample"

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    path '*_fastqc.zip', emit: zip
    path '*_fastqc.html', emit: html

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${reads[0]} ${reads[1]}
    """
}


process ALIGN_READS {

    tag "$sample"

    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path indexed_reference

    output:
    tuple val(sample), path("${sample}.bam"),
        emit: bam

    script:
    """
    set -o pipefail

    bwa mem \
        -t ${task.cpus} \
        -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' \
        ${indexed_reference}/reference.fa \
        ${reads[0]} ${reads[1]} \
    | samtools sort \
        -o ${sample}.bam -

    samtools quickcheck -v ${sample}.bam
    """
}


process ALIGNMENT_QC {

    tag "$sample"

    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    path "${sample}.bam.bai", emit: index
    path "${sample}.flagstat.txt", emit: flagstat
    path "${sample}.stats.txt", emit: stats
    path "${sample}.coverage.txt", emit: coverage

    script:
    """
    samtools index ${bam} ${sample}.bam.bai

    samtools flagstat ${bam} \
        > ${sample}.flagstat.txt

    samtools stats ${bam} \
        > ${sample}.stats.txt

    samtools coverage ${bam} \
        > ${sample}.coverage.txt
    """
}


process MULTIQC {

    tag 'Combined sequencing QC'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports

    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    script:
    """
    multiqc . --force
    """
}


workflow {

    samples = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample,
                [
                    file(row.r1, checkIfExists: true),
                    file(row.r2, checkIfExists: true)
                ]
            )
        }

    reference = channel.value(
        file(
            params.reference,
            checkIfExists: true
        )
    )

    INDEX_REFERENCE(reference)

    FASTQC(samples)

    ALIGN_READS(
        samples,
        INDEX_REFERENCE.out
    )

    ALIGNMENT_QC(
        ALIGN_READS.out.bam
    )

    reports = FASTQC.out.zip
        .mix(
            ALIGNMENT_QC.out.flagstat,
            ALIGNMENT_QC.out.stats,
            ALIGNMENT_QC.out.coverage
        )
        .collect()

    MULTIQC(reports)
}
