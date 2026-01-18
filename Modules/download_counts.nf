process DOWNLOAD_COUNTS {

  tag { sample_id }
  publishDir "${params.outdir}/01_raw_counts", mode: 'copy'

  input:
  tuple val(sample_id), val(counts)

  output:
  tuple val(sample_id), path("${sample_id}.csv")

  script:
  """
  set -euo pipefail

  if [[ "${counts}" =~ ^https?:// ]]; then
    echo "[INFO] Downloading: ${counts}"
    curl -L --retry 5 --retry-delay 5 -o ${sample_id}.csv.gz "${counts}" || wget -O ${sample_id}.csv.gz "${counts}"
    gzip -dc ${sample_id}.csv.gz > ${sample_id}.csv
  else
    echo "[INFO] Using local file: ${counts}"
    if [[ "${counts}" == *.gz ]]; then
      gzip -dc "${counts}" > ${sample_id}.csv
    else
      cp "${counts}" ${sample_id}.csv
    fi
  fi
  """
}
