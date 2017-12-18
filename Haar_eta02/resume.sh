#!/bin/bash

#matlab -nodisplay -nodesktop -r "run init_calc.m"

anz=$(<anz.txt)
j=$(<j.txt)

while (( j <= anz ))
    do
	echo "bash: j =" $j

	echo " ----------------------------------------------- "
	echo " ----------------- Call Matlab ----------------- "
	echo " ----------------------------------------------- "

	matlab -nodisplay -nodesktop -r "run calc.m"

	echo " ----------------------------------------------- "
	echo " ---------------- Close Matlab ----------------- "
	echo " ----------------------------------------------- "
	
	j=$(<j.txt)
done


