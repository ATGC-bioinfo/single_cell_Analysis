process SCANPY_DE_AND_PLOTS {
  tag { sample_id }
  publishDir "${params.outdir}/02_analysis/${sample_id}/06_scanpy_de_and_plots", mode: 'copy'

  input:
  tuple val(sample_id), path(in_h5ad), path(model_dir), path(figures_in), path(tables_in)

  output:
  tuple val(sample_id),
        path("figures"),
        path("tables")

  script:
  """
  python3 ${projectDir}/data/06_scanpy_de_and_plots.py \
    --sample_id "${sample_id}" \
    --in_h5ad "${in_h5ad}" \
    --figdir figures \
    --tabledir tables
  """
}
