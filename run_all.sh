#!/usr/bin/env bash
# Convenience wrapper: full WCT-motivated GWTC log-domain diagnostic pipeline.
# Diagnostic only -- does not prove WCT, does not replace LVK population inference.
set -euo pipefail
cd "$(dirname "$0")"

NULL_PRIMARY="${NULL_PRIMARY:-1000}"
NULL_STRESS="${NULL_STRESS:-400}"
NK="${NK:-120}"

echo "== fetch =="          ; python3 scripts/fetch_gwtc_catalog.py --offline-ok
echo "== build =="          ; python3 scripts/build_gwtc_table.py
echo "== primary scan =="   ; python3 scripts/run_gwtc_log_scan.py --variable M_chirp --cumulative --p-astro-min 0.5 --bins 20 --null-n "$NULL_PRIMARY" --n-k "$NK" --seed 12345
echo "== nulls dump =="     ; python3 scripts/run_gwtc_nulls.py --variable M_chirp --cumulative --p-astro-min 0.5 --bins 20 --null-n "$NULL_PRIMARY" --n-k "$NK" --seed 12345
echo "== bin stress =="     ; python3 scripts/run_gwtc_bin_stress.py --variable M_chirp --cumulative --p-astro-min 0.5 --null-n "$NULL_STRESS" --n-k "$NK"
echo "== variable stress ==" ; python3 scripts/run_gwtc_variable_stress.py --cumulative --p-astro-min 0.5 --bins 20 --null-n "$NULL_STRESS" --n-k "$NK" --include-bounded
echo "== threshold stress ==" ; python3 scripts/run_gwtc_threshold_stress.py --variable M_chirp --bins 20 --null-n "$NULL_STRESS" --n-k "$NK"
echo "== controls =="       ; python3 scripts/run_gwtc_controls.py --variable M_chirp --cumulative --p-astro-min 0.5 --bins 20 --null-n "$NULL_STRESS" --n-k "$NK"
echo "== master =="         ; python3 scripts/make_gwtc_master_table.py
echo "== verdict =="        ; python3 scripts/make_gwtc_verdict.py --primary-variable M_chirp --primary-subset cumulative_pa0.5
echo "== done -- see outputs/summary/VERDICT.txt =="
