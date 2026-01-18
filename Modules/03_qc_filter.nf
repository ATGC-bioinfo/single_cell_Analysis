process QC_FILTER {
  tag { sample_id }
  publishDir "${params.outdir}/02_analysis/${sample_id}/03_qc_filter", mode: 'copy'

  input:
  tuple val(sample_id), path(in_h5ad)

  output:
  tuple val(sample_id),
        path("${sample_id}.qc.h5ad"),
        path("figures"),
        path("tables")

  script:
  """
  python3 ${projectDir}/data/03_qc_filter.py \
    --sample_id "${sample_id}" \
    --in_h5ad "${in_h5ad}" \
    --out_h5ad "${sample_id}.qc.h5ad" \
    --figdir figures \
    --tabledir tables
  """
}
