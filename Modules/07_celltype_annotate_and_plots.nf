process CELLTYPE_ANNOTATE_AND_PLOTS {
  tag { sample_id }
  publishDir "${params.outdir}/02_analysis/${sample_id}/07_celltype_annotate_and_plots", mode: 'copy'

  input:
  tuple val(sample_id), path(in_h5ad), path(fig_de), path(tab_de)

  output:
  path "${sample_id}.final.h5ad", emit: final_h5ad
  path "figures", emit: figures_dir
  path "tables", emit: tables_dir

  script:
  """
  python3 ${projectDir}/data/07_celltype_annotate_and_plots.py \
    --sample_id "${sample_id}" \
    --in_h5ad "${in_h5ad}" \
    --out_h5ad "${sample_id}.final.h5ad" \
    --figdir figures \
    --tabledir tables
  """
}
