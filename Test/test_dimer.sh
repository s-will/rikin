#!/bin/bash

RRIKDIR=/home/will/Research/Projects/RRIkinetics/RRIKin
TDIR=$RRIKDIR/Test
TFILE=$TDIR/dimer

COMMON_OPTS=""
COMMON_OPTS="$COMMON_OPTS" #  --max-hyb-length 8
ENUM_OPTS="--max-hyb-length-diff 5 --max-hyb-energy 6 --max-total-energy 10"

#seqA=CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAG
seqA=CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAG
#seqA=CGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGAGCUAGCUGAAACGUAGCAGACAGAUACUAGUUUCAAAACUUCCUUGACAGAACCAUUUCAUGUUCAAUAAUGAAAAUUACUUUCACAUGUUUUAGUGGAAAACGUACGUACGUAUCGUAGCGGUUGGACUUACGUAUAC

#seqB=GCUGCGAUGCAUGCCUCGAU
seqB=GCUGCGAUGCAUGCCUCGAUAUGCCUCGAUAUGCCUCGAU
#seqB=GCUGCGAUGCAUGCCUCGAUAUGCCUCGAUAUGCCUCGAUGUAAAAAUUAAAAAAAUAUAAUACAUUAAAUGCAAAAUAAGUUUAGCUUACGGUACGUAG


## nette Beispiele:
# competing sites
#seqA=AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA
#seqB=CCCGCCC


function call {
    echo
    echo $*
    time $*
}

function call_redirect {
    tgt=$1
    shift 1
    echo
    echo $* \> $tgt
    time $* > $tgt
}

call_redirect $TFILE $RRIKDIR/src/rrikin_enum $seqA $seqB --no-double --verbose $ENUM_OPTS $COMMON_OPTS

sort -g -k1 $TFILE > $TFILE.s

call $RRIKDIR/src/rrikin_barriers $seqA $seqB <$TFILE.s -o $TFILE.bg --no-double --verbose --max-recover-energy 0   --min-rate 1e-9 $COMMON_OPTS

call_redirect $TFILE.out $RRIKDIR/src/rrikin_prune $TFILE.bg --pffile $TFILE.pf  --min-rate 1e-9 --qu-out 10 --qu-pequ 75 --num-out 1000 --num-pequ 1000 --verbose

call src/xrates.m $TFILE.pf --out $TFILE.kin --t8 1e15

call Rscript $RRIKDIR/src/kinetics.R $TFILE.kin $TDIR/kinetics.pdf
