#!/bin/bash

# rikin_pipeline.sh
# Run the Rikin pipeline

# -----
# USAGE
# -----
usage() {
    cat <<+++USAGE_MESSAGE
============================================================
rikin_pipeline.sh

Run the Rikin pipeline


USAGE: rikin_pipeline.sh [-h|--help] [--jobid jobname] [--config file] SEQA SEQB

EXAMPLE: rikin_pipeline.sh -j test CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGACAGAUACUAGUUUCAAAACUUCCUUGACAGAACCAUUUCAUGUUCAAUAAUGAAAAUUACUUUCACAUGUUUUAGUGGAAAACGUACGUACGUAUCGUAGCGGUUGGACUUACGUAUAC GCUGCGAUGCAUGCCUCGAUAUGCCUCGAUAUGCCUCGAUGUAAAAAUUAAAAAAAUAUAAUACAUUAAAUGCAAAAUAAGUUUAGCUUACGGUACGUAG

--jobid defines identifier for the job; it must be a valid file and directory name; output is written to files jobname/jobname*

--config If the given config file or rikin_pipeline.cfg exists in the current directory, then the it defines the pipeline configuration. It must be a valid bash file that defines COMMON_OPTS and ENUM_OPTS. Attention: this script is executed.

SEQA, SEQB: RNA sequences as words over A,C,G,U
============================================================

+++USAGE_MESSAGE
}


# -----------------
# Default arguments
# -----------------
#

# jobname and name of output directory
JOBID="Rikin.results"
CONFIGURATION=""

COMMON_OPTS="--max-hyb-length-diff 8 --span=100"
ENUM_OPTS="--max-hyb-energy 6 --max-total-energy 10"

# ---------------
# Parse arguments
# ---------------

VALID_ARGS=$(getopt -o hj:c: --long help,jobid:,configuration: -- "$@")
if [[ $? -ne 0 ]]; then
    echo "Argument parsing failed."
    usage
    exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -h | --help)
        usage
        exit 0
        ;;
    -j | --jobid)
        JOBID="$2"
        shift 2
        ;;
    -c | --configuration)
        CONFIGURATION="$2"
        shift 2
        ;;
    --)
        shift;
        break
        ;;

    :)
      echo "Option requires an argument."
      usage
      exit 1
      ;;

    ?)
      echo -e "Invalid command option."
      usage
      exit 1
      ;;
  esac
done


SEQA=$1
SEQB=$2
shift 2

if [ "$SEQA" == "" -o "$SEQB" == "" ]; then
    echo "No sequences specified"
    usage
    exit -1
fi


# ------------------
# Check arg validity
# ------------------
#
# existence of files ...


if [ -e "$JOBID" ] ; then
    echo "WARNING: output directory $JOBID exists already."
else
    mkdir "$JOBID"
fi
FILENAME="$JOBID/$JOBID"


if [ "$CONFIGURATION" == "" ]; then
    for CNAME in "rikin_pipeline.cfg"; do
        if [ -e $CNAME ]; then
            CONFIGURATION="$CNAME"
        fi
    done
fi

if [ "$CONFIGUATION" != "" ]; then
    if [ -e "$CONFIGURATION" ]; then
        . "$CONFIGURATION"
        cp "$CONFIGURATION" "$JOBID"
    else
        echo "ERROR: configuration file $CONFIGURATION does not exist."
        exit -1
    fi
fi

# ----------
# Print input summary
# ----------


cat <<+++INPUT_SUMMARY
============================================================
RIkin Pipeline

Job ID:      $JOBID
SeqA:        $SEQA
SeqB:        $SEQB
COMMON_OPTS: $COMMON_OPTS
ENUM_OPTS:   $ENUM_OPTS

Start date:  $(date)
Working dir: $(pwd)
============================================================
+++INPUT_SUMMARY


# ----------------
# MAIN

function call {
    echo
    echo $*
    $*
}

function tcall {
    echo
    echo $*
    time $*
}

function call_redirect {
    tgt=$1
    shift 1
    echo
    echo $* \> $tgt
    $* > $tgt
}


call_redirect $FILENAME rikin_enum $SEQA $SEQB --no-double --verbose $ENUM_OPTS $COMMON_OPTS

LC_ALL=C sort -g -k1 $FILENAME > $FILENAME.s

call rikin_barriers $SEQA $SEQB -i $FILENAME.s -o $FILENAME.bg --track $FILENAME.barriers_track.gz --no-double --verbose $COMMON_OPTS

# call_redirect $FILENAME.out $RIKDIR/_inst/bin/rikin_prune $FILENAME.bg --pffile $FILENAME.pf  --min-rate 1e-10 --num-out 50 --num-pequ 100 --track $FILENAME.prune_track.gz --verbose --barfile $FILENAME.bar --ratesfile $FILENAME.rates

call_redirect $FILENAME.out rikin_prune $FILENAME.bg --pffile $FILENAME.pf  --min-rate 1e-7 --num-out 10 --num-pequ 60 --track $FILENAME.prune_track.gz --verbose --barfile $FILENAME.bar --ratesfile $FILENAME.rates


tcall rikin_xrates.m $FILENAME.pf --out $FILENAME.kin --t8 1e15 --mfpts

call rikin_kinetics.R $FILENAME.kin -o $FILENAME.pdf


cat <<+++END_MESSAGE


============================================================
RIkin pipeline finished at $(date).

Output files written to directory $JOBID
============================================================
+++END_MESSAGE
