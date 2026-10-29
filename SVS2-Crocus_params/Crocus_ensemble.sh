#!/bin/bash
set -uo pipefail
# NOTE: -e is intentionally NOT set globally anymore. With an ensemble loop,
# we want a single crashed member to be logged and skipped, not to kill the
# whole run. Error checking for critical commands is done explicitly below.

# Define option sets
HSNOWDRIFT=(R21W R21R R21F)
LSNOWDRIFT=(.true.)
HSNOWMETAMO=(B21)
HSNOWFALL=(R21 GW1 GW2)
HSNOWCOND=(S97)
HSNOWRES=(DEF)
# Ensemble sets
HSNOWHOLD=(B92 B02 O04 SPK)
HSNOWCOMP=(B92 S14 T11) 

# CHANGE THE RECORDS FILE NAME EVERY TIME YOU RUN THE SCRIPT !!!IMPORTANT!!!
sim_name=CROCUS_ENSEMBLE_NAME

RecordFile=log_$sim_name.text

# --- Directory layout ---------------------------------------------------
# work_dir      : ephemeral scratch space used *while a member is running*.
#                 Wiped and rebuilt every iteration, exactly as before, so
#                 sa_mesh always starts from a clean directory.
# output_link   : exp_dir/output -> work_dir. This is the symlink sa_mesh /
#                 generate_nc_output.py expect to find during/after a run.
#                 Behaviour unchanged from the original script.
#
# collection_dir: NEW - a single, PERSISTENT folder (lives in the work area,
#                 never wiped) that collects every renamed output .nc file
#                 from every ensemble member.
# collection_link: NEW - exp_dir/all_outputs -> collection_dir, so from the
#                 experiment folder you have one shortcut to browse every
#                 result, instead of the old scattered per-combo folders.
# --------------------------------------------------------------------------

work_dir=$sim_name
exp_dir=$sim_name
output_link="$exp_dir/output"

collection_dir=${sim_name}_outputs
collection_link="$exp_dir/all_outputs"

log_dir="$exp_dir/logs"
crash_dir="$log_dir/crashed"

mkdir -p "$exp_dir" "$log_dir" "$crash_dir" "$collection_dir"

# Persistent shortcut to the collection folder - set up once, not per loop.
ln -sfn "$collection_dir" "$collection_link"

# Initialise counter and write header
c=1
echo "c HSNOWDRIFT LSNOWDRIFT HSNOWMETAMO HSNOWFALL HSNOWCOND HSNOWHOLD HSNOWCOMP HSNOWRES" > "$RecordFile"

# Loop through ALL combinations
for drift in "${HSNOWDRIFT[@]}"; do
  for ldrift in "${LSNOWDRIFT[@]}"; do
    for metam in "${HSNOWMETAMO[@]}"; do
      for fall in "${HSNOWFALL[@]}"; do
        for cond in "${HSNOWCOND[@]}"; do
          for hold in "${HSNOWHOLD[@]}"; do
            for comp in "${HSNOWCOMP[@]}"; do
              for res in "${HSNOWRES[@]}"; do

                echo "Loop iteration: $c"
                echo "$c $drift $ldrift $metam $fall $cond $hold $comp $res" >> "$RecordFile"

                # NOTE: this naming pattern only encodes hold/comp/fall/drift.
                # If you ever give ldrift/metam/cond/res more than one value
                # each, two different combinations could produce the same
                # name and the "already exists" check below would wrongly
                # skip a genuinely new run. Add the missing variables to the
                # name if/when those arrays grow.
                expnamefolder="CaSR32R_norain_HOLD${hold}_COMP${comp}_FALL${fall}_DRIFT${drift}"
                final_nc="$collection_dir/${expnamefolder}.nc"

                echo "$c $drift $ldrift $metam $fall $cond $hold $comp $res"

                # --- Skip already-completed members --------------------
                if [[ -f "$final_nc" ]]; then
                  echo "  Output already exists ($final_nc) - skipping."
                  ((c++))
                  continue
                fi

                echo "    Executing model"

                sed -e "s/ens_drift/$drift/g" \
                    -e "s/ens_ldrift/$ldrift/g" \
                    -e "s/ens_metamo/$metam/g" \
                    -e "s/ens_fall/$fall/g" \
                    -e "s/ens_cond/$cond/g" \
                    -e "s/ens_hold/$hold/g" \
                    -e "s/ens_comp/$comp/g" \
                    -e "s/ens_res/$res/g" \
                    MESH_parameters_orig.txt > MESH_parameters.txt

                # --- Rebuild the work dir and the output symlink from scratch ---
                # This guarantees "output" can never be a leftover real directory,
                # which is what makes ln -sf nest a symlink inside it instead of
                # replacing it (the actual cause of the recursive folders).
                rm -rf "$work_dir"
                mkdir -p "$work_dir"
                ln -sfn "$work_dir" "$output_link"

                echo " MESH parameters changed, ready to run sa_mesh"

                # --- Fail-safe: don't let one crashed member kill the ensemble ---
                sim_log="$log_dir/output_${c}.log"
                sa_mesh_ok=1
                if ! ./sa_mesh > "$sim_log" 2>&1; then
                  sa_mesh_ok=0
                fi

                if [[ $sa_mesh_ok -eq 0 ]]; then
                  crash_log="$crash_dir/CRASHED_${expnamefolder}.log"
                  cp -f "$sim_log" "$crash_log"
                  echo "  !! sa_mesh crashed for $expnamefolder (exit code $?) - see $crash_log"
                  ((c++))
                  continue
                fi

                echo "sa_mesh complete"

                # nc generation can also fail even if sa_mesh itself didn't
                nc_ok=1
                if ! python /code/generate_nc_output.py svs2 $sim_name; then
                  nc_ok=0
                fi

                if [[ $nc_ok -eq 0 || ! -f "$output_link/out_svs2.nc" ]]; then
                  crash_log="$crash_dir/CRASHED_${expnamefolder}.log"
                  cp -f "$sim_log" "$crash_log"
                  echo "  !! nc generation failed/missing for $expnamefolder - see $crash_log"
                  ((c++))
                  continue
                fi

                # --- Move the renamed output into the single collection folder ---
                cp "$output_link/out_svs2.nc" "$final_nc"
                cp -f "$RecordFile" "/work/$RecordFile"

                ((c++))
              done
            done
          done
        done
      done
    done
  done
done