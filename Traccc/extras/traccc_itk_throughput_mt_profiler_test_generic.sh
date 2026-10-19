#!/bin/bash
#
# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#
# Generic script running the selected instance of the multi-threaded throughput
# executable on a whole set of ITk ttbar simulations, with different pileup
# values. Supports multiple architectures and data types.
#

usage() {
   echo "Script running a suite of multi-threaded throughput tests"
   echo ""
   echo "Usage: traccc_throughput_mt_profiler_test_generic.sh [options]"
   echo ""
   echo "Basic options:"
   echo "  -a <architecture>    Architecture (ada or other) [default: ada]"
   echo "  -T <type>            Data type (fs, jets, muon) [default: fs]"
   echo "  -x <executable>      Selects the executable to use"
   echo "  -i <inputDir>        Selects the input directory with all files"
   echo "  -m <minThreads>      Minimum number of threads to test"
   echo "  -t <maxThreads>      Maximum number of threads to test"
   echo "  -s <threadStep>      Steps to increase the thread count by"
   echo "  -r <repetitions>     The number of repetitions in the test"
   echo "  -e <eventMultiplier> Multiplier for the number of events per thread"
   echo "  -c <csvFile>         Name of the output CSV file"
   echo "  -y <throughputType>  Type of throughput test to run (traccc/g200/g100)"
   echo "  -h                   Print this help"
   echo ""
}

# Parse command line arguments
while getopts ":a:T:x:i:m:t:s:r:c:y:h" opt; do
   case $opt in
      a)
         ARCHITECTURE=$OPTARG
         ;;
      T)
         TYPE=$OPTARG
         ;;
      x)
         TRACCC_EXECUTABLE=$OPTARG
         ;;
      i)
         TRACCC_INPUT_DIR=$OPTARG
         ;;
      m)
         TRACCC_MIN_THREADS=$OPTARG
         ;;
      t)
         TRACCC_MAX_THREADS=$OPTARG
         ;;
      s)
         TRACCC_THREAD_STEP=$OPTARG
         ;;
      r)
         TRACCC_REPETITIONS=$OPTARG
         ;;
      c)
         TRACCC_CSV_FILE=$OPTARG
         ;;
      y)
         TRACCC_THROUGPUT_TYPE=$OPTARG
         ;;
      h)
         usage
         exit 0
         ;;
      :)
         echo "Argument -$OPTARG requires a parameter!"
         usage
         exit 1
         ;;
      ?)
         echo "Unknown argument: -$OPTARG"
         usage
         exit 1
         ;;
   esac
done

# Determine executable path based on architecture
if [[ "${ARCHITECTURE:-ada}" == "ada" ]]; then
   DEFAULT_EXECUTABLE="/afs/cern.ch/user/e/exochell/traccc/acts_ada/acts/traccc-build/bin/traccc_throughput_mt_cuda"
else
   DEFAULT_EXECUTABLE="/afs/cern.ch/user/e/exochell/traccc/acts/acts/traccc-build/bin/traccc_throughput_mt_cuda"
fi

# Set default values for all configuration options
TRACCC_ARCHITECTURE=${ARCHITECTURE:-"ada"}
TRACCC_TYPE=${TYPE:-"fs"}
TRACCC_EXECUTABLE=${TRACCC_EXECUTABLE:-"${DEFAULT_EXECUTABLE}"}
TRACCC_INPUT_DIR=${TRACCC_INPUT_DIR:-"/eos/project/a/atlas-eftracking/GPU/ITk_data/traccc_standalone_data"}
TRACCC_MIN_THREADS=${TRACCC_MIN_THREADS:-1}
TRACCC_MAX_THREADS=${TRACCC_MAX_THREADS:-12}
TRACCC_THREAD_STEP=${TRACCC_THREAD_STEP:-1}
TRACCC_REPETITIONS=${TRACCC_REPETITIONS:-5}
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
TRACCC_OUTPUT_DIR=${TRACCC_OUTPUT_DIR:-"/eos/user/e/exochell/traccc/profile_output/${TRACCC_TYPE}/${TRACCC_ARCHITECTURE}_gbts"}
TRACCC_CSV_FILE=${TRACCC_CSV_FILE:-"${TRACCC_TYPE}_${TRACCC_ARCHITECTURE}_output_${TIMESTAMP}.csv"}
TRACCC_THROUGPUT_TYPE=${TRACCC_THROUGPUT_TYPE:-"traccc"}

# Determine input data directory based on type
case "${TRACCC_TYPE}" in
   fs)
      TRACCC_DATA_DIR="${TRACCC_INPUT_DIR}/ttbar_mu200"
      ;;
   jets)
      TRACCC_DATA_DIR=${TRACCC_DATA_DIR:-"/eos/user/e/exochell/traccc/traccc_athena_plots/g200/traccc-athena/data/roiInputJets"}
      ;;
   muon)
      TRACCC_DATA_DIR=${TRACCC_DATA_DIR:-"/eos/user/e/exochell/traccc/traccc_athena_plots/g200/traccc-athena/data/roiInputMuon"}
      ;;
   *)
      echo "***"
      echo "*** Unknown type: '${TRACCC_TYPE}'"
      echo "*** Supported types: fs, jets, muon"
      echo "***"
      exit 1
      ;;
esac

# Print the configuration received.
echo "Using configuration:"
echo "   ARCHITECTURE    : ${TRACCC_ARCHITECTURE}"
echo "   TYPE            : ${TRACCC_TYPE}"
echo "   EXECUTABLE      : ${TRACCC_EXECUTABLE}"
echo "   INPUT_DIR       : ${TRACCC_INPUT_DIR}"
echo "   DATA_DIR        : ${TRACCC_DATA_DIR}"
echo "   MIN_THREADS     : ${TRACCC_MIN_THREADS}"
echo "   MAX_THREADS     : ${TRACCC_MAX_THREADS}"
echo "   THREAD_STEP     : ${TRACCC_THREAD_STEP}"
echo "   REPETITIONS     : ${TRACCC_REPETITIONS}"
echo "   CSV_FILE        : ${TRACCC_OUTPUT_DIR}/${TRACCC_CSV_FILE}"
echo "   THROUGHPUT_TYPE : ${TRACCC_THROUGPUT_TYPE}"

# Check whether the output file already exists. Refuse to overwrite existing files.
if [[ -f "${TRACCC_OUTPUT_DIR}/${TRACCC_CSV_FILE}" ]]; then
   echo "***"
   echo "*** Will not overwrite ${TRACCC_CSV_FILE}!"
   echo "***"
   exit 1
fi

# Additional flags tuning the cuts for the G100/G200 pipelines.
G200_CUTS=(--seedfinder-z-range=-3000.:3000.
           --seedfinder-r-range=33.:320.
           --seedfinder-vertex-range=-114.:114.
           --seedfinder-minPt=0.9
           --seedfinder-cotThetaMax=27.2899
           --seedfinder-deltaR-range=20.:100.
           --seedfinder-deltaZMax=800.
           --seedfinder-impactMax=2.
           --seedfinder-sigmaScattering=3.
           --seedfinder-maxPtScattering=10.
           --seedfinder-maxSeedsPerSpM=2
           --max-num-branches-per-seed=3
           --max-num-branches-per-surface=1
           --track-candidates-range=7:20
           --min-step-length-for-next-surface=0.5
           --max-step-counts-for-next-surface=100
           --chi2-max=30.
           --max-num-skipping-per-cand=2
           --stepping-min-stepsize=0.0001
           --rk-tolerance-mm=0.0001
           --stepping-path-limit=5.
           --stepping-max-rk-updates=10000
           --stepping-use-mean-loss=1
           --stepping-use-eloss-gradient=0
           --stepping-use-field-gradient=0
           --stepping-do-covariance-transport=1
           --overstep-tolerance-um=-300.
           --max-num-tracks-per-measurement 1
           --initial-links-per-seed 20
           --seedfinder-compatSeedLimit 1)

G100_CUTS=${G200_CUTS[@]}
G100_CUTS+=(--reco-stage=seeding)

# Select which of these flags to use.
TRACCC_CUTS=()
if [[ "${TRACCC_THROUGPUT_TYPE}" == "g200" ]]; then
   TRACCC_CUTS=${G200_CUTS[@]}
elif [[ "${TRACCC_THROUGPUT_TYPE}" == "g100" ]]; then
   TRACCC_CUTS=${G100_CUTS[@]}
elif [[ "${TRACCC_THROUGPUT_TYPE}" != "traccc" ]]; then
   echo "***"
   echo "*** Unknown throughput type: '${TRACCC_THROUGPUT_TYPE}'"
   echo "***"
   exit 1
fi

# Check the full path
if [[ -f "${TRACCC_OUTPUT_DIR}/${TRACCC_CSV_FILE}" ]]; then
   echo "***"
   echo "*** Will not overwrite ${TRACCC_OUTPUT_DIR}/${TRACCC_CSV_FILE}!"
   echo "***"
   exit 1
fi

# Put a header on the CSV file.
echo "directory,threads,loaded_events,cold_run_events,processed_events,warm_up_time,processing_time" \
   > "${TRACCC_OUTPUT_DIR}/${TRACCC_CSV_FILE}"

# Counter for a nice printout.
COUNTER=1
COUNT=$((${TRACCC_MAX_THREADS}*${TRACCC_REPETITIONS}))

# Iterate over the number of threads.
for NTHREAD in $(seq ${TRACCC_MIN_THREADS} ${TRACCC_THREAD_STEP} ${TRACCC_MAX_THREADS}); do
   # Perform the requested number of repetitions.
   for REPEAT in $(seq ${TRACCC_REPETITIONS}); do
      echo ""
      echo "Running test ${COUNTER} / ${COUNT}"
      ((COUNTER++))

      ${TRACCC_EXECUTABLE}                                                  \
         --detector-file="${TRACCC_INPUT_DIR}/ITk_DetectorBuilder_geometry.json" \
         --material-file="${TRACCC_INPUT_DIR}/ITk_detector_material.json"   \
         --grid-file="${TRACCC_INPUT_DIR}/ITk_DetectorBuilder_surface_grids.json" \
         --digitization-file="${TRACCC_INPUT_DIR}/ITk_digitization_config.json" \
         --conditions-file="${TRACCC_INPUT_DIR}/ITk_conditions_config.json" \
         --read-bfield-from-file                                            \
         --bfield-file="${TRACCC_INPUT_DIR}/ITk_bfield.cvf"                 \
         --input-directory="${TRACCC_DATA_DIR}"                             \
         --use-acts-geom-source=0                                           \
         --input-events=10                                                  \
         --cpu-threads=${NTHREAD}                                           \
         --cold-run-events=$((5*${NTHREAD}))                                \
         --processed-events=$((100*${NTHREAD}))                             \
         --log-file="${TRACCC_OUTPUT_DIR}/${TRACCC_CSV_FILE}"               \
         --input-skip=5                                                     \
         --useGBTS \
         --gbts_config_dir="${TRACCC_INPUT_DIR}" \
         ${TRACCC_CUTS[@]}
   done
done
