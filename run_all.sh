#!/usr/bin/env bash
# Runs the four checks in the order used for expected_output.txt.
set -e
cd "$(dirname "$0")/scripts"
for f in check_background_and_horndeski_map.py check_distance_ratio_numerical.py check_tensor_action.py fit_standard_templates.py; do
  echo "== $f"
  python3 "$f"
done
