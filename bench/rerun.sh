#!/usr/bin/env bash
# Re-render pipeline stages in order, retrying the ones that depend on a remote repository.
#   usage: bash bench/rerun.sh [stage ...]
#   stages: massive metabolights workbench 02 03 04 05   (default: 02 03 04 05)
#
# Sequential by design: BiocFileCache is single-writer, so two renders at once deadlock.
# GNPS2 and the EBI FTP drop requests under load, which looks like a hard failure but is
# not, so a stage that fails on a known-transient message is retried rather than abandoned.
# Anything else stops the run with the error printed, because a silent retry loop around a
# real bug wastes hours.
cd "$(dirname "$0")/.." || exit 1

export CONFIRM_CAP="${CONFIRM_CAP:-Inf}"   # confirmation is uncapped unless told otherwise
TRANSIENT='Failed to connect to GNPS2|Gateway Timeout|Failed to perform HTTP|bfcrpath|database is locked|Connection reset|Timeout was reached|Failed to connect to MetaboLights'

declare -A QMD=(
  [massive]=confirmation/massIVE-hit-confirmation
  [metabolights]=confirmation/metaboLights-confirmation
  [workbench]=confirmation/metabolomics-workbench-confirmation
  [02]=02-confirmation-and-bio-files
  [03]=03-all-datasets
  [04]=04-cyp-validation
  [05]=05-meta-analysis
)
declare -A TRIES=( [massive]=40 [metabolights]=20 [workbench]=20 [02]=5 [03]=3 [04]=3 [05]=3 )

stages=("${@:-02 03 04 05}"); stages=(${stages[@]})

for st in "${stages[@]}"; do
  q="${QMD[$st]}"
  [ -n "$q" ] || { echo "unknown stage: $st"; exit 2; }
  log="bench/rerun-$(basename "$q").log"
  max="${TRIES[$st]:-3}"

  for attempt in $(seq 1 "$max"); do
    echo "STAGE $(date +%H:%M) $st attempt $attempt/$max"
    t0=$(date +%s)
    quarto render "$q.qmd" > "$log" 2>&1
    rc=$?
    echo "STAGE $(date +%H:%M) $st exit=$rc after $(( ($(date +%s) - t0) / 60 )) min"
    [ $rc -eq 0 ] && break

    if grep -qiE "$TRANSIENT" "$log" 2>/dev/null; then
      echo "STAGE $(date +%H:%M) transient failure, waiting 10 min"
      sleep 600
    else
      echo "STAGE $(date +%H:%M) $st failed:"
      grep -B2 -A10 -iE "^Error|Error in" "$log" | head -30
      exit 1
    fi
  done
  [ $rc -eq 0 ] || { echo "STAGE $st exhausted $max attempts"; exit 1; }
done
echo "STAGE ALL COMPLETE $(date +%H:%M)"
