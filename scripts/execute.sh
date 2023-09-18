#!/bin/bash
DIR_ALGO="/home/user/Documents/algorithms"
DIR_D="/home/user/Documents/dataset/benchmark/labelled/email-Enron"
DIR_Q="/home/user/Documents/dataset/benchmark/queries/query5"

START=$(date +%s.%N)

#Executing RapidMatch
${DIR_ALGO}/RapidMatch/build/matching/RapidMatch.out -d ${DIR_D}RM.csv -q ${DIR_Q}RM.csv -order nd -num MAX #> /dev/null 2>&1

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
printf "############# Computed version: ${DIFF}s\n\n"

START=$(date +%s.%N)
#Executing CECI
${DIR_ALGO}/ceci-release/ceci ${DIR_D}CECI.csv ${DIR_Q}CECI.csv #> /dev/null 2>&1
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
printf "############# Computed version: ${DIFF}s\n\n"

START=$(date +%s.%N)
#Executing DAF
${DIR_ALGO}/DAF/bin/daf_parallel_10min -d ${DIR_D}DAF.csv -q ${DIR_Q}DAF.csv -n 1 -m 1000000000 -h 8 #> /dev/null 2>&1
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
printf "############# Computed version: ${DIFF}s\n\n"



