nextflow.enable.dsl=2
params.bam       = null
params.reference = null
params.truth     = null
params.bed       = null
params.region    = 'chr20:10000000-11000000'
params.outdir    = 'runs'

process REGIONAL_CALL {
  tag 'HG002-chr20'
  publishDir params.outdir, mode: 'copy'
  input:
  path bam
  path ref
  path fai
  path dictfile
  path truth
  path truth_index
  path bed
  val region
  output:
  path 'report.html'
  path 'benchmark.json'
  path 'variants.norm.vcf.gz'
  script:
  """
  mkdir -p data results
  ln -s ${bam} data/sample.bam
  ln -s ${ref} data/chr20.fa
  ln -s ${fai} data/chr20.fa.fai
  ln -s ${dictfile} data/chr20.dict
  ln -s ${truth} data/truth.vcf.gz
  ln -s ${truth_index} data/truth.vcf.gz.tbi
  ln -s ${bed} data/truth.bed
  cp -r ${projectDir}/scripts ./scripts
  REGION='${region}' BAM=data/sample.bam REF=data/chr20.fa TRUTH=data/truth.vcf.gz BED=data/truth.bed bash scripts/run.sh
  cp results/report.html results/benchmark.json results/variants.norm.vcf.gz .
  """
}

workflow {
  if (!params.bam || !params.reference || !params.truth || !params.bed)
    error 'Specify --bam --reference --truth --bed'
  REGIONAL_CALL(file(params.bam), file(params.reference), file("${params.reference}.fai"),
    file(params.reference.replaceFirst(/\.[^.]+$/, '.dict')),
    file(params.truth), file("${params.truth}.tbi"), file(params.bed), params.region)
}
