process SOLO_DOUBLET {
  tag { sample_id }
  publishDir "${params.outdir}/02_analysis/${sample_id}/02_doublets", mode: 'copy'

  input:
  tuple val(sample_id), path(in_h5ad)

  output:
  tuple val(sample_id),
        path("${sample_id}.nodoublet.h5ad"),
        path("${sample_id}.doublets.csv")

  script:
  """
  python3 ${projectDir}/data/02_doublet_solo.py \
    --in_h5ad "${in_h5ad}" \
    --out_h5ad "${sample_id}.nodoublet.h5ad" \
    --out_doublets "${sample_id}.doublets.csv" \
    --max_epochs ${params.max_epochs ?: 400}
  """
}
