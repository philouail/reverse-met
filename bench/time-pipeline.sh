#!/usr/bin/env bash
# Time each analysis-notebook render, one at a time (no contention), and write a
# throughput table. Extraction is a separate one-off step; its per-dataset
# wall-clock is logged to peaks_review_v2/_extract_times.csv by bench/peak_review_panels.R.
#   usage: bash bench/time-pipeline.sh [01 02 03 04 05]
cd "$(dirname "$0")/.." || exit 1
OUT=artifacts/pipeline_timing.csv
stages=("${@:-01 02 03 04 05}"); stages=(${stages[@]})
echo "stage,notebook,seconds,exit" > "$OUT"
declare -A QMD=( [01]=01-masst-curation [02]=02-confirmation-and-bio-files
                 [03]=03-all-datasets [04]=04-meta-analysis [05]=05-cyp-validation )
for st in "${stages[@]}"; do
  q="${QMD[$st]}"; [ -f "$q.qmd" ] || continue
  extra=""; [ "$st" = "03" ] && extra="-P render_reports:true"
  t0=$(date +%s)
  quarto render "$q.qmd" $extra > "bench/render_$st.log" 2>&1; rc=$?
  t1=$(date +%s)
  echo "$st,$q,$((t1-t0)),$rc" | tee -a "$OUT"
done
echo "wrote $OUT"
