#!/usr/bin/env bash
# bash bench/rerun.sh <application> [stage ...]: re-renders stages in order, retrying network
# failures. Sequential: BiocFileCache is single-writer, so two renders at once deadlock.
ROOT="$(cd "$(dirname "$0")/.." && pwd)" || exit 1
APP="${1:?usage: bash bench/rerun.sh <application> [stage ...]}"; shift
cd "$ROOT/application/$APP" 2>/dev/null || { echo "no such application: $APP"; exit 2; }

export CONFIRM_CAP="${CONFIRM_CAP:-Inf}"   # confirmation is uncapped unless told otherwise
# Notebooks stop on an unreachable repository so this script retries them.
# NETWORK must match NETWORK_ERROR in R/metadata.R word for word.
NETWORK='Failed to perform HTTP request|resolve host|connect to server|Timeout (of [0-9]+ seconds )?was reached|timed out|Connection (was )?reset|Connection refused|Recv failure|Send failure|receiving data|sending data|Server returned nothing|Transferred a partial file|SSL connect error|(HTTP|error:?) ?(429|50[234])|Too Many Requests|Bad Gateway|Service Unavailable|Gateway Time-?out|not all .rnames. found'
TRANSIENT="Failed to connect to GNPS2|Failed to connect to MetaboLights|bfcrpath|database is locked|$NETWORK"

declare -A QMD=(
  [massive]=confirmation/massIVE-hit-confirmation
  [metabolights]=confirmation/metaboLights-confirmation
  [workbench]=confirmation/metabolomics-workbench-confirmation
  [01]=01-search
  [02]=02-metadata
  [03]=03-curation
  [04]=04-confirmation
  [05]=05-extraction
  [06]=06-analysis
)
declare -A TRIES=( [massive]=40 [metabolights]=20 [workbench]=20 [01]=3 [02]=5 [03]=3 [04]=5 [05]=3 [06]=3 )

stages=("${@:-04 05 06}"); stages=(${stages[@]})

for st in "${stages[@]}"; do
  q="${QMD[$st]}"
  [ -n "$q" ] || { echo "unknown stage: $st"; exit 2; }
  [ -f "$q.qmd" ] || { echo "STAGE $st skipped: $APP has no $q.qmd"; continue; }
  log="$ROOT/bench/rerun-$APP-$(basename "$q").log"
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
