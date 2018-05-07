#!/bin/bash
experimentDir=$PWD/experiment

ARCH=(avx2 "avx2 sr" avx512f "avx512f sml")
ALPHA_VALUES=("0.1" "0.1,0.75,1.5")
CATEGORIES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
PINVAR=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
PMATRIX_ITR=(10000 20000 30000 40000 50000 60000 70000 80000 90000 100000)
SITES=(10000 20000 30000 40000 50000 60000 70000 80000 90000 100000)
REPEAT=(1 2 3 4 5 6 7 8 9 10)

mkdir -p $experimentDir

for arch in "${ARCH[@]}"
do
    echo "Benchmark for Exec Mode: $arch"
    for alphaValues in "${ALPHA_VALUES[@]}"
    do
        subAlphaValues=${alphaValues//,/__}
        pinvarDir=$experimentDir/pinvar
        mkdir -p ${pinvarDir}
        for pinvar in "${PINVAR[@]}"
        do
            outputDir=${pinvarDir}/${arch/ /_}__${pinvar/./_}__${subAlphaValues//./_}
            mkdir -p $outputDir
            for repeat in "${REPEAT[@]}"
            do
                ./obj/derivatives-benchmark ${arch} -alpha=${alphaValues} -p-invar=${pinvar} > $outputDir/${repeat}.txt
            done
        done

        sitesDir=$experimentDir/sites
        mkdir -p ${sitesDir}
        for sites in "${SITES[@]}"
        do
            outputDir=${sitesDir}/${arch/ /_}__${sites}__${subAlphaValues//./_}
            mkdir -p $outputDir
            for repeat in "${REPEAT[@]}"
            do
                ./obj/derivatives-benchmark ${arch} -alpha=${alphaValues} -n-sites=${sites} > $outputDir/${repeat}.txt
            done
        done

        categoriesDir=$experimentDir/categories
        mkdir -p ${categoriesDir}
        for categories in "${CATEGORIES[@]}"
        do
            outputDir=${categoriesDir}/${arch/ /_}__${categories}__${subAlphaValues//./_}
            mkdir -p $outputDir
            for repeat in "${REPEAT[@]}"
            do
                ./obj/derivatives-benchmark ${arch} -alpha=${alphaValues} -n-categories=${categories} > $outputDir/${repeat}.txt
            done
        done
    done

    pmatrixDir=$experimentDir/pmatrix
    mkdir -p ${pmatrixDir}
    for pmatrix in "${PMATRIX_ITR[@]}"
    do
        outputDir=${pmatrixDir}/${arch/ /_}__${pmatrix}
        mkdir -p $outputDir
        for repeat in "${REPEAT[@]}"
        do
            ./obj/derivatives-benchmark ${arch} -n-pmatrix-itr=${pmatrix} -n-sites=1000000 > $outputDir/${repeat}.txt
        done
    done
done

echo "Done"
