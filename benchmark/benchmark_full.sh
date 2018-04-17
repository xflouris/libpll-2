#!/bin/bash
experimentDir=$PWD/experiment

ARCH=(avx2 "avx2 sr" avx512f "avx512f sml")
ALPHA_VALUES=("0.1" "0.1,0.75,1.5")
CATEGORIES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
PINVAR=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.9)
PMATRIX_ITR=(1 100 200 300 400 500 600 700 800 900 1000)
SITES=(1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000)

mkdir -p $experimentDir

echo "Full Benchmark"
pinvarDir=$experimentDir/full
mkdir -p ${pinvarDir}
for arch in "${ARCH[@]}"
do
    for pinvar in "${PINVAR[@]}"
    do
        for alphaValues in "${ALPHA_VALUES[@]}"
        do
            for categories in "${CATEGORIES[@]}"
            do
                for sites in "${SITES[@]}"
                do
                    subAlphaValues=${alphaValues//,/__}
                    fileName=${pinvarDir}/${arch/ /_}__${pinvar/./_}__${subAlphaValues//./_}__${categories}__${sites}.txt
                    ./obj/derivatives-aa-benchmark ${arch} -alpha=${alphaValues} -p-invar=${pinvar} -n-categories=${categories} -n-sites=${sites} > $fileName
                done
            done
        done
    done
done
