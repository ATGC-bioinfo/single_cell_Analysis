#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { DOWNLOAD_COUNTS }               from './Modules/download_counts.nf'
include { MAKE_ADATA }                    from './Modules/01_make_adata.nf'
include { SOLO_DOUBLET }                  from './Modules/02_doublet_solo.nf'
include { QC_FILTER }                     from './Modules/03_qc_filter.nf'
include { SCANPY_CLUSTER }                from './Modules/04_scanpy_cluster.nf'
include { SCVI_INTEGRATE }                from './Modules/05_scvi_integrate.nf'
include { SCANPY_DE_AND_PLOTS }           from './Modules/06_scanpy_de_and_plots.nf'
include { CELLTYPE_ANNOTATE_AND_PLOTS }   from './Modules/07_celltype_annotate_and_plots.nf'

workflow {

  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample_id as String, row.counts as String) }
    .set { samples_ch }

  downloaded = DOWNLOAD_COUNTS(samples_ch)

  adata0 = MAKE_ADATA(downloaded)

  adata1 = SOLO_DOUBLET(adata0)

  adata2 = QC_FILTER(adata1.map { sample_id, h5ad_path, doublets_csv -> tuple(sample_id, h5ad_path) })

  adata3 = SCANPY_CLUSTER(adata2)

  adata4 = SCVI_INTEGRATE(adata3)

  de_out = SCANPY_DE_AND_PLOTS(adata4)

  // Extract h5ad from SCVI_INTEGRATE and join with DE results
  scvi_h5ad = adata4.map { sample_id, h5ad_path, scvi_model_path, figures_scvi, tables_scvi ->
    tuple(sample_id, h5ad_path)
  }
  
  // Join h5ad with DE results: scvi_h5ad is (sample_id, h5ad_path), de_out is (sample_id, figures, tables)
  joined_for_final = scvi_h5ad.join(de_out, by: 0)
  
  CELLTYPE_ANNOTATE_AND_PLOTS(
    joined_for_final.map { joined_list ->
      // join returns: [sample_id, h5ad_path, figures, tables]
      tuple(joined_list[0], joined_list[1], joined_list[2], joined_list[3])
    }
  )
}
