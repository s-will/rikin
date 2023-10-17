#!/bin/bash

# rikin_pipeline.sh
# Run the Rikin pipeline

## version of this script -- will be included in ouptut
VERSION="0.9"


# -----
# USAGE
# -----
usage() {
    cat <<+++USAGE_MESSAGE
============================================================
rikin_pipeline.sh ver $VERSION

Run the Rikin pipeline


USAGE: rikin_pipeline.sh [-h|--help] [options] SEQA SEQB

EXAMPLE CALL:
rikin_pipeline.sh -j example AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA CCCGCCC 2>&1 | tee example.out

--outdir output directory
 
--config The given config file defines the pipeline configuration. It must be a valid bash file that is sourced such that it can define environment variables to control the pipeline. These comprise, COMMON_OPTS, ENUM_OPTS. The file is sourced after the global configuration file rikin_pipeline.cfg.

--dryrun don't run commands and/or write files

SEQA, SEQB: RNA sequences as words over A,C,G,U
============================================================

+++USAGE_MESSAGE
}

# -----------------
# Bin dir
BINDIR=$(which $0)
BINDIR=$(dirname $BINDIR)

# -----------------
# Default arguments
# -----------------
# jobname and name of output directory
OUTDIR="RIkin-pipeline-results"
GLOBAL_CONFIGURATION="rikin_pipeline.cfg"
CONFIGURATION=""
DRYRUN=false # off
REUSE=false #off

SEQA_NAME=seqA
SEQB_NAME=seqB

# ---------------
# Parse arguments
# ---------------

VALID_ARGS=$(getopt -o hj:c:o: --long help,jobid:,configuration:,outdir:,dryrun,reuse -- "$@")
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
    -o | --outdir)
        OUTDIR="$2"
        shift 2
        ;;
    -j | --jobid)
        JOBID="$2"
        shift 2
        ;;
    -c | --configuration)
        CONFIGURATION="$2"
        shift 2
        ;;
    --dryrun)
        DRYRUN=true
        shift
        ;;

    --reuse)
        REUSE=true
        shift
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

if [ ! -e "$OUTDIR" ] ; then
    mkdir "$OUTDIR"
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
# Source configuration
# ----------

# a) global config
GLOBAL_CONFIGURATION_FULL="$BINDIR/$GLOBAL_CONFIGURATION"
if [ -e "$GLOBAL_CONFIGURATION_FULL" ] ; then
    . "$GLOBAL_CONFIGURATION_FULL"
else
    echo "WARNING: Cannot find global config at \"$GLOBAL_CONFIGURATION_FULL\"."
fi

# b) optional config
if [ "$CONFIGUATION" != "" ]; then
    . "$CONFIGURATION"
fi


# ----------
# Print input summary
# ----------


cat <<+++INPUT_SUMMARY
============================================================
RIkin Pipeline ver $VERSION

Input:
  Arguments: $*
  Job ID:    $JOBID
  SeqA:      $SEQA
  SeqB:      $SEQB
  Outdir:    $OUTDIR

Environment:
  Working dir: $(pwd)
  BINDIR=$BINDIR
  RIkin version $($BINDIR/rikin_enum --version)

Options:
  COMMON_OPTS="$COMMON_OPTS" (controls rikin_enum and rikin_barriers)
  ENUM_OPTS="$ENUM_OPTS" (controls rikin_enum)
  BARRIERS_OPTS="$BARRIERS_OPTS" (controls rikin_barriers)
  PRUNE_OPTS="$PRUNE_OPTS" (controls rikin_prune)
  XRATES_OPTS="$XRATES_OPTS" (conrols solving of master equation by rikin_xrates.m)
  KINETICS_OPTS="$KINETICS_OPTS" (controls plotting by rikin_kinetics.R)

Start date:  $(date)
============================================================
+++INPUT_SUMMARY


# ----------------
# MAIN

function call {
    echo
    echo $*
    if ! $DRYRUN ; then
        $*
    fi
}

function tcall {
    echo
    echo $*
    if ! $DRYRUN ; then
        time $*
    fi
}

function call_redirect {
    tgt=$1
    shift 1
    echo
    echo $* \> $tgt
    if ! $DRYRUN ; then
        $* > $tgt
    fi
}

########################################
## enumerate and sort
if [ -e ${OUTDIR}/sorted_states ] \
   ; then
    \gzip -f ${OUTDIR}/sorted_states
fi

if ! $REUSE || [ ! -e ${OUTDIR}/sorted_states.gz ] \
   ; then

    call_redirect ${OUTDIR}/states $BINDIR/rikin_enum \
		  $SEQA $SEQB \
		  $COMMON_OPTS $ENUM_OPTS

    if [ $? -ne 0 ] ; then
	echo ERROR: rikin_enum died
	rm ${OUTDIR}/sorted_states.gz >/dev/null 2>&1
	exit -1
    fi

    LC_ALL=C sort -k1,1n ${OUTDIR}/states > ${OUTDIR}/sorted_states
    rm ${OUTDIR}/states

    REUSE=false
fi

########################################
## barriers
if ! $REUSE \
       || [ ! -e ${OUTDIR}/bg ] \
       || [ ! -e ${OUTDIR}/barriers_track.gz ] \
       || [ ! -e ${OUTDIR}/track-ipps-barriers.gz ] \
   ; then

    if [ -e ${OUTDIR}/sorted_states.gz ] \
       ; then
	\gunzip -f ${OUTDIR}/sorted_states.gz
    fi

    call $BINDIR/rikin_barriers $SEQA $SEQB \
	 -i ${OUTDIR}/sorted_states \
	 -o ${OUTDIR}/bg \
	 --compress-track --track ${OUTDIR}/barriers_track.gz \
	 --track-ipps=${OUTDIR}/track-ipps-barriers.gz \
	 $COMMON_OPTS $BARRIERS_OPTS
    if [ $? -ne 0 ] ; then
	echo ERROR: rikin_barriers died
	rm ${OUTDIR}/bg >/dev/null 2>&1
	exit -1
    fi
    REUSE=false
fi

if [ -e ${OUTDIR}/sorted_states ] \
   ; then
    \gzip -f ${OUTDIR}/sorted_states
fi

########################################
## prune
if ! $REUSE \
       || [ ! -e ${OUTDIR}/out ] \
       || [ ! -e ${OUTDIR}/pf ] \
       || [ ! -e ${OUTDIR}/bar ] \
       || [ ! -e ${OUTDIR}/rates.gz ] \
       || [ ! -e ${OUTDIR}/track-ipps-prune.gz  ] \
   ; then
    call_redirect ${OUTDIR}/out $BINDIR/rikin_prune ${OUTDIR}/bg \
		  $PRUNE_OPTS \
		  --pffile ${OUTDIR}/pf --verbose --barfile ${OUTDIR}/bar \
                  --ratesfile ${OUTDIR}/rates.gz \
		  --compress-track --track ${OUTDIR}/prune_track.gz \
		  --track-pps-out=${OUTDIR}/track-ipps-prune.gz \
		  --track-pps-in=${OUTDIR}/track-ipps-barriers.gz
    if [ $? -ne 0 ] ; then
	echo ERROR: rikin_prune died
	rm ${OUTDIR}/rates.gz >/dev/null 2>&1
	exit -1
    fi
    REUSE=false
fi

########################################
## solve master equation
if ! $REUSE \
       || [ ! -e ${OUTDIR}/kin ] \
   ; then
    tcall $BINDIR/rikin_xrates.m ${OUTDIR}/pf --out ${OUTDIR}/kin $XRATES_OPTS
    if [ $? -ne 0 ] ; then
	echo >${OUTDIR}/error "$THE_CALL"
	echo ERROR: xrates died
	rm ${OUTDIR}/kin >/dev/null 2>&1
	exit -1
    else
	touch ${OUTDIR}/xrates.done
	rm ${OUTDIR}/error >/dev/null 2>&1
    fi
    REUSE=false
fi

########################################
## generate plot
if ! $REUSE \
       || [ ! -e ${OUTDIR}/kinetics.pdf ] \
   ; then
    call $BINDIR/rikin_kinetics.R ${OUTDIR}/kin -o ${OUTDIR}/kinetics.pdf \
	 --pps ${OUTDIR}/track-ipps-prune.gz \
	 --seqA $SEQA --seqB $SEQB --nameA ${SEQA_NAME} --nameB ${SEQB_NAME} \
         $KINETICS_OPTS
    if [ $? -ne 0 ] ; then
	echo ERROR: kinetics.R died
	rm ${OUTDIR}/kinetics.pdf
	exit -1
    fi
    REUSE=false
fi


## -------------------------------------------
# 
cat <<+++END_MESSAGE


============================================================
RIkin pipeline finished at $(date).

Output files written to directory $JOBID
============================================================
+++END_MESSAGE
