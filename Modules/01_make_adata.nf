process MAKE_ADATA {
  tag { sample_id }
  publishDir "${params.outdir}/02_analysis/${sample_id}/01_make_adata", mode: 'copy'

  input:
  tuple val(sample_id), path(counts_csv)

  output:
  tuple val(sample_id), path("${sample_id}.raw.h5ad")

  script:
  """
  python3 ${projectDir}/data/01_make_adata.py \
    --sample_id "${sample_id}" \
    --counts_csv "${counts_csv}" \
    --out_h5ad "${sample_id}.raw.h5ad"
  """
}
