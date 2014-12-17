#!/bin/bash

RRIKDIR=/home/will/Research/Projects/RRIkinetics/RRIKin
TDIR=$RRIKDIR/Test
TFILE=$TDIR/dimer

COMMON_OPTS=""
COMMON_OPTS="$COMMON_OPTS --max-hyb-length-diff 5" #  --max-hyb-length 8
ENUM_OPTS="--max-hyb-energy 6 --max-total-energy 10"

seqA=CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAG
#seqA=CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAG
#seqA=CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGACAGAUACUAGUUUCAAAACUUCCUUGACAGAACCAUUUCAUGUUCAAUAAUGAAAAUUACUUUCACAUGUUUUAGUGGAAAACGUACGUACGUAUCGUAGCGGUUGGACUUACGUAUAC

seqB=GCUGCGAUGCAUGCCUCGAU
#seqB=GCUGCGAUGCAUGCCUCGAUAUGCCUCGAUAUGCCUCGAU
#seqB=GCUGCGAUGCAUGCCUCGAUAUGCCUCGAUAUGCCUCGAUGUAAAAAUUAAAAAAAUAUAAUACAUUAAAUGCAAAAUAAGUUUAGCUUACGGUACGUAG


## nette Beispiele:
# competing sites
#seqA=AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA
#seqB=CCCGCCC

## starker decoy state
#seqA=CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGACAGAUACUAGUUUCAAAACUUCCUUGACAGAACCAUUUCAUGUUCAAUAAUGAAAAUUACUUUCACAUGUUUUAGUGGAAAACGUACGUACGUAUCGUAGCGGUUGGACUUACGUAUAC
#seqB=GCUGCGAUGCAUGCCUCGAU


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

call_redirect $TFILE $RRIKDIR/src/rrikin_enum $seqA $seqB --no-double --verbose $ENUM_OPTS $COMMON_OPTS

sort -g -k1 $TFILE > $TFILE.s

call $RRIKDIR/src/rrikin_barriers $seqA $seqB -i $TFILE.s -o $TFILE.bg --track $TFILE.barriers_track --no-double --verbose $COMMON_OPTS

# call_redirect $TFILE.out $RRIKDIR/src/rrikin_prune $TFILE.bg --pffile $TFILE.pf  --min-rate 1e-10 --num-out 50 --num-pequ 100 --track $TFILE.prune_track --verbose --barfile $TFILE.bar --ratesfile $TFILE.rates

call_redirect $TFILE.out $RRIKDIR/src/rrikin_prune $TFILE.bg --pffile $TFILE.pf  --min-rate 1e-7 --num-out 10 --num-pequ 60 --track $TFILE.prune_track --verbose --barfile $TFILE.bar --ratesfile $TFILE.rates


tcall src/xrates.m $TFILE.pf --out $TFILE.kin --t8 1e15 --mfpts

call Rscript $RRIKDIR/src/kinetics.R $TFILE.kin $TDIR/kinetics.pdf
