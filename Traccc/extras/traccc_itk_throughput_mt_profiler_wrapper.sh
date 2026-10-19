#!/bin/bash
#
# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#
# Master wrapper script for multi-threaded throughput tests.
# Calls the generic profiler script with specified architecture and type.
#


#example:   bash ./traccc_itk_throughput_mt_profiler_wrapper.sh -a ada -T muon -m 1 -t 3 -r 2
#           bash ./traccc_itk_throughput_mt_profiler_wrapper.sh -a ada -T jets -m 1 -t 1 -r 10
#           bash ./traccc_itk_throughput_mt_profiler_wrapper.sh -a tesla -T fs -m 1 -t 1 -r 10

usage() {
   echo "Master wrapper for TRACCC throughput tests"
   echo ""
   echo "Usage: traccc_itk_throughput_mt_profiler_test_wrapper.sh -a <arch> -T <type> [options]"
   echo ""
   echo "Required options:"
   echo "  -a <architecture>    Architecture (ada or other)"
   echo "  -T <type>            Data type (fs, jets, muon)"
   echo ""
   echo "Optional options (passed to generic script):"
   echo "  -x <executable>      Selects the executable to use"
   echo "  -i <inputDir>        Selects the input directory with all files"
   echo "  -m <minThreads>      Minimum number of threads to test"
   echo "  -t <maxThreads>      Maximum number of threads to test"
   echo "  -s <threadStep>      Steps to increase the thread count by"
   echo "  -r <repetitions>     The number of repetitions in the test"
   echo "  -c <csvFile>         Name of the output CSV file"
   echo "  -y <throughputType>  Type of throughput test to run (traccc/g200/g100)"
   echo "  -h                   Print this help"
   echo ""
}

ARCHITECTURE=""
TYPE=""

# Parse architecture and type first
while getopts ":a:T:h" opt; do
   case $opt in
      a)
         ARCHITECTURE=$OPTARG
         ;;
      T)
         TYPE=$OPTARG
         ;;
      h)
         usage
         exit 0
         ;;
      *)
         ;;
   esac
done

if [[ -z "${ARCHITECTURE}" ]] || [[ -z "${TYPE}" ]]; then
   echo "Error: Both -a (architecture) and -T (type) are required!"
   echo ""
   usage
   exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENERIC_SCRIPT="${SCRIPT_DIR}/traccc_itk_throughput_mt_profiler_test_generic.sh"

if [[ ! -f "${GENERIC_SCRIPT}" ]]; then
   echo "Error: Generic profiler script not found at ${GENERIC_SCRIPT}"
   exit 1
fi

# Call the generic script with architecture and type, passing all original arguments
exec "${GENERIC_SCRIPT}" -a "${ARCHITECTURE}" -T "${TYPE}" "$@"
