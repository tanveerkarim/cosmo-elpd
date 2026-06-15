#!/bin/bash

FILE=$1

if [ -z "$FILE" ]; then
    echo "Usage: $0 <cobaya_output_file>"
    exit 1
fi

echo "===================================================="
echo "           COBAYA / POLYCHORD DIAGNOSTIC"
echo "===================================================="
echo "File: $FILE"
echo "----------------------------------------------------"

# --- core counts ---
logl=$(grep -c "Computed log-likelihood" "$FILE")
camb_ok=$(grep -c "Got input parameters" "$FILE")
planck=$(grep -c "planck_2018" "$FILE")

camb_fail=$(grep -c "HMCode INTEGRATE\|Computation of cosmological products failed" "$FILE")

posterior=$(grep -c "Posterior to be computed" "$FILE")
live=$(grep -c "nlive\|live points" "$FILE")

lock=$(grep -c -i "locked" "$FILE")

# --- CAMB failure fraction ---
if [ $((camb_ok + camb_fail)) -gt 0 ]; then
    camb_frac=$(awk -v f=$camb_fail -v o=$camb_ok 'BEGIN {printf "%.5f", f/(f+o)}')
else
    camb_frac="NA"
fi

# --- ratios ---
if [ "$planck" -gt 0 ]; then
    camb_planck_ratio=$(awk -v c=$camb_ok -v p=$planck 'BEGIN {printf "%.3f", c/p}')
else
    camb_planck_ratio="NA"
fi

# --- print table ---
printf "\n%-30s %10s\n" "METRIC" "VALUE"
printf "%-30s %10s\n" "--------------------------------" "----------"

printf "%-30s %10d\n" "LogL evaluations" "$logl"
printf "%-30s %10d\n" "CAMB evaluations" "$camb_ok"
printf "%-30s %10d\n" "Planck likelihood hits" "$planck"
printf "%-30s %10d\n" "Posterior evaluations" "$posterior"
printf "%-30s %10d\n" "Live-point events" "$live"
printf "%-30s %10d\n" "Lock-file hits" "$lock"
printf "%-30s %10d\n" "CAMB failures" "$camb_fail"

printf "%-30s %10s\n" "CAMB failure fraction" "$camb_frac"
printf "%-30s %10s\n" "CAMB/Planck ratio" "$camb_planck_ratio"

echo "----------------------------------------------------"

# --- quick interpretation hints ---
echo "INTERPRETATION:"
echo "- CAMB failure fraction > 0.01  => HMCode stress regime"
echo "- CAMB/Planck ratio << 1        => Planck dominates cost"
echo "- Lock-file hits > 0            => potential MPI filesystem contention"
echo "- High Planck counts            => CamSpec dominates runtime"
echo "----------------------------------------------------"
