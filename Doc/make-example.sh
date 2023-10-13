#!/bin/bash

rikin_pipeline.sh -j example AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA CCCGCCC | tee example.out

convert example/example.pdf example.png

if [ ! -e html/Doc ] ; then mkdir html/Doc ; fi
cp example.png html/Doc

