process SCANPY_CLUSTER {
  tag { sample_id }
  publishDir "${params.outdir}/02_analysis/${sample_id}/04_scanpy_cluster", mode: 'copy'

  input:
  tuple val(sample_id), path(in_h5ad), path(figures_in), path(tables_in)

  output:
  tuple val(sample_id),
        path("${sample_id}.scanpy.h5ad"),
        path("figures"),
        path("tables")

  script:
  """
  python3 ${projectDir}/data/04_scanpy_cluster.py \
    --sample_id "${sample_id}" \
    --in_h5ad "${in_h5ad}" \
    --out_h5ad "${sample_id}.scanpy.h5ad" \
    --figdir figures \
    --tabledir tables
  """
}
