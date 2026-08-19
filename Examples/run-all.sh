#!/usr/bin/env sh

# RNAInterKin: run all examples

# Simple shell script to run the RNAInterKin pipeline for kinetcics
# calculations and generation of plots for all of the RNA--RNA interactions
# examples in the Examples directory.
#
# USAGE: bash run-all.sh
#
# Run this script from the Examples directory.
#
# The run produces output on directories example, HP1-HP2, HP3-HP2, MicA-MicA, and MicA-OmpA.
# These directories contain results and intermediary data of the kinetics
# computation as well as plots for visual analysis.
# (see Documentation for the output of rikin_pipeline.py)
#

echo "Run rikin_pipeline.py on the RNAInterKin examples in the Examples directory"

# --------------------
# Toy example
rikin_pipeline.py -o example --seqA AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA --seqB CCCGCCC

# --------------------
# HP1:HP2
rikin_pipeline.py --fastaA HP1.fasta --fastaB HP2.fasta -c khp.cfg -o HP1-HP2

# --------------------
# HP3:HP2
rikin_pipeline.py --fastaA HP3.fasta --fastaB HP2.fasta -c khp.cfg -o HP3-HP2


# --------------------
# MicA-MicA homo-dimer
rikin_pipeline.py --fastaA MicA.fasta --fastaB MicA.fasta -c MicA-MicA.cfg -o MicA-MicA

# --------------------
# MicA-OmpA hetero-dimer
rikin_pipeline.py --fastaA MicA.fasta --fastaB OmpA.fasta -c MicA-ompA.cfg -o MicA-OmpA
