#!/bin/bash

rikin_pipeline.sh -o example AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA CCCGCCC 2>&1 | tee example.out

convert example/kinetics.pdf example.png

if [ ! -e html/Doc ] ; then mkdir html/Doc ; fi
cp example.png html/Doc

