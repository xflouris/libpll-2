#!/bin/bash
experimentDir=$PWD/experiment

ARCH=(avx2 "avx2 sr" avx512f "avx512f sml")
ALPHA_VALUES=("0.1" "0.1,0.75,1.5")
CATEGORIES=(1 5 25)
PINVAR=(0.0 0.3 0.5 0.6)
SITES=(1000 1000 100000 1000000)

mkdir -p $experimentDir

otherDir=$experimentDir/other
mkdir -p ${otherDir}

for arch in "${ARCH[@]}"
do
    for pinvar in "${PINVAR[@]}"
    do
        #subAlphaValues=${alphaValues//,/__}
        #fileName=${arch/ /_}__${pinvar/./_}__${subAlphaValues//./_}__${categories}__${sites}
        fileName=${otherDir}/${arch/ /_}__${pinvar/./_}.txt
        ./obj/derivatives-aa-benchmark ${arch} -p-invar=${pinvar} > $fileName
    done
done

#mkdir -p $MAIN_DIR/detail
#./obj/derivatives-aa-benchmark-simple avx2        > $MAIN_DIR/detail/avx2.txt
#./obj/derivatives-aa-benchmark-simple avx2 sr     > $MAIN_DIR/detail/avx2_sr.txt
#./obj/derivatives-aa-benchmark-simple avx512f     > $MAIN_DIR/detail/avx512.txt
#./obj/derivatives-aa-benchmark-simple avx512f sml > $MAIN_DIR/detail/avx512_sml.txt
#action.py $MAIN_DIR/detail/*.txt     > $MAIN_DIR/detail/detail_summary.csv
