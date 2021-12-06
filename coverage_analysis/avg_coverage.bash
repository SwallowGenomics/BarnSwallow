#!/bin/bash

while read line
do
awk -v scaffold="${line}"  '{if($1==scaffold) print $0}' $2 > ${line}.genomecov

if ! [[ -e ${line}.genomecov_500bp ]] ; then

awk -v last="$(wc -l ${line}.genomecov)" '{sum+=$3; counter+=1; if ((NR%500==0) || (NR==last-1)){print $1"\t"$2-counter"\t"$2"\t"sum/counter; sum=0; counter=0}}' ${line}.genomecov > ${line}.genomecov_500bp

else

echo "${line}.genomecov_500bp already exists"

fi

done<$1