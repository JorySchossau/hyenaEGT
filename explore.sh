#!/bin/bash

TRIAL=trial
EXECUTABLE=public
REPLICATES=4
REP_START=1

R_LOW=3
R_HIGH=7

N_THREADS=8
ZETA=5

EXEC_STATUS=""
mkdir ${TRIAL}
cd ${TRIAL}
cp ../${EXECUTABLE} .
counter=1
for each_rep in $(seq ${REP_START} 1 ${REPLICATES}); do
	for each_r in $(seq ${R_LOW} 1 ${R_HIGH}); do
		echo "starting simulation for ${each_rep}.${each_r}"
		if [[ ${counter} == ${N_THREADS} || ( ${each_rep} == ${REPLICATES} && ${each_r} == ${R_HIGH} ) ]]; then
			./${EXECUTABLE} --experiment trial --replicate ${each_rep} --end out.${each_rep}_${each_r}.end --beta 0 --gamma 0 --fcdmi 0.25 0.25 0.25 0.25 --updates 150000 --zeta ${ZETA} --r ${each_r}
			counter=0
		else
			./${EXECUTABLE} --experiment trial --replicate ${each_rep} --end out.${each_rep}_${each_r}.end --beta 0 --gamma 0 --fcdmi 0.25 0.25 0.25 0.25 --updates 150000 --zeta ${ZETA} --r ${each_r} &
		fi
		counter=$((counter+1))
	done
done
sleep 3
