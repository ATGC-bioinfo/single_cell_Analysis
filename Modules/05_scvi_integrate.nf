process SCVI_INTEGRATE {
  tag { sample_id }
  publishDir "${params.outdir}/02_analysis/${sample_id}/05_scvi_integrate", mode: 'copy'

  input:
  tuple val(sample_id), path(in_h5ad), path(figures_in), path(tables_in)

  output:
  tuple val(sample_id),
        path("${sample_id}.scvi.h5ad"),
        path("scvi_model"),
        path("figures"),
        path("tables")

  script:
  """
  python3 ${projectDir}/data/05_scvi_integrate.py \
    --sample_id "${sample_id}" \
    --in_h5ad "${in_h5ad}" \
    --out_h5ad "${sample_id}.scvi.h5ad" \
    --model_dir scvi_model \
    --figdir figures \
    --tabledir tables \
    --max_epochs ${params.max_epochs ?: 400}
  """
}
