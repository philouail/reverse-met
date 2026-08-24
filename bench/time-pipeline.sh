#!/usr/bin/env bash
# Time each analysis-notebook render, one at a time.
#   usage: bash bench/time-pipeline.sh [01 03 04 05] [--reports]
# --reports also regenerates the per-dataset peak-review sheets in stage 03.
# Writes artifacts/pipeline_timing.csv.
cd "$(dirname "$0")/.." || exit 1

REPORTS=false
args=()
for a in "$@"; do
  if [ "$a" = "--reports" ]; then REPORTS=true; else args+=("$a"); fi
done
set -- "${args[@]}"

OUT=artifacts/pipeline_timing.csv
stages=("${@:-01 03 04 05}"); stages=(${stages[@]})
declare -A QMD=( [01]=01-masst-curation [02]=02-confirmation-and-bio-files
                 [03]=03-all-datasets [04]=04-cyp-validation
                 [05]=05-meta-analysis )

VERS=$(Rscript -e 'cat(paste0("R ",R.version$major,".",R.version$minor,", Chromatograms ",utils::packageVersion("Chromatograms"),", Spectra ",utils::packageVersion("Spectra"),", xcms ",utils::packageVersion("xcms")))' 2>/dev/null)
{
  echo "# run: $(date +%Y-%m-%d) on $(nproc) cores, $(uname -s)"
  echo "# $VERS"
  echo "# quarto $(quarto --version 2>/dev/null)"
  echo "# reports: $REPORTS"
  echo "stage,notebook,seconds,exit"
} > "$OUT"

for st in "${stages[@]}"; do
  q="${QMD[$st]}"; [ -f "$q.qmd" ] || { echo "skip $st (no $q.qmd)"; continue; }
  extra=""; [ "$st" = "03" ] && [ "$REPORTS" = true ] && extra="-P render_reports:true"
  echo "timing $st ($q) ..."
  t0=$(date +%s)
  quarto render "$q.qmd" $extra > "bench/render_$st.log" 2>&1; rc=$?
  t1=$(date +%s)
  printf '%s,%s,%s,%s\n' "$st" "$q" "$((t1-t0))" "$rc" | tee -a "$OUT"
done
echo "wrote $OUT"
