#!/bin/bash

d1=$(date +"%s")

rm -rf output
rm -rf log
mkdir -p output
mkdir -p log
inputDir=$(ls ./input_pa1)
echo ${inputDir}
for i_dat in ${inputDir}
do
	echo "---------------------------------------------------------"
	echo "Now executing file:  "${i_dat}
	./bin/fm ./input_pa1/${i_dat} ./output/${i_dat}.txt > ./log/${i_dat}.log
	echo ${i_dat}" finished!!"
	echo "---------------------------------------------------------"
done

d2=$(date +"%s")
echo "total used $((${d2}-${d1})) seconds"

exit 0
