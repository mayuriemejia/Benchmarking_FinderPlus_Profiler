process AMRFINDER_DB_UPDATE {
  /*
    # (•˕ •マ.ᐟ
    ✧. *✧. *✧. *✧.
    AMRFINDER_DB_UPDATE · snapshot de base de datos
    ✧. *✧. *✧. *✧.
  */

  tag { "amrfinder-db" }
  publishDir { "${params.outdir}/04_amrfinder/db" }, mode: 'copy'
  container params.img_amr

  output:
  path('amrfinder_db')

  script:
  """
  set -euo pipefail
  mkdir -p amrfinder_db

  if [[ -n "${params.amrfinder_db}" && -d "${params.amrfinder_db}" ]]; then
      echo "[INFO] Usando base de datos local: ${params.amrfinder_db}"
      cp -r ${params.amrfinder_db}/* amrfinder_db/
  else
      echo "[INFO] Descargando base de datos AMRFinderPlus..."
      update_out=\$(amrfinder -u 2>&1 || true)
      dbdir=\$(echo "\$update_out" | awk -F': ' '/^Database directory:/{print \$2}' | tail -n 1)
      [[ -z "\$dbdir" ]] && dbdir="/usr/local/bin/data"
      [[ ! -d "\$dbdir" ]] && { echo "DB no encontrada: \$dbdir" >&2; exit 1; }
      cp -r "\$dbdir/"* amrfinder_db/
  fi
  """
}
