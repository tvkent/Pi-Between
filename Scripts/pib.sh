#!/bin/bash

#python pib.py

cd ../Results

files=(ai si)
for file in files
do
	for n in (2..12)
	do
	sed -i '1d' ${file}_${n}.txt
	done
done


cat si_1.txt si_2.txt si_3.txt si_4.txt si_5.txt si_6.txt si_7.txt si_8.txt si_9.txt si_10.txt si_11.txt si_12.txt > si.txt
cat ai_1.txt ai_2.txt ai_3.txt ai_4.txt ai_5.txt ai_6.txt ai_7.txt ai_8.txt ai_9.txt ai_10.txt ai_11.txt ai_12.txt > ai.txt

