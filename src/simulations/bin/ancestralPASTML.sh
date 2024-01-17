#!/bin/bash
params=$1
matrix=$2
msalen=$3
msatsv=$4
tree=$5
threads=$6
outprefix=$7

# Sort list of parameter files by site number
VAR=$(tr ' ' '\n' <<<"$params" | sort -V | paste -sd' ' -)
# Keep only the site indices
ID=$(echo $VAR | sed 's/parameter_//g' | tr ' ' '\n' | sort -V | paste -sd' ')
# Print the matrix name msalen times (to give to pastml)
#MATRIX=$(printf "$matrix %.0s" {1..$msalen})
MATRIX=$(yes "$matrix" | head -n $msalen | tr '\n' ' ')

gunzip -c ${msatsv} > pastml.input

# Call pastml
#echo pastml --threads $threads -t ${tree} -d pastml.input --prediction_method MAP -m CUSTOM_RATES --rate_matrix $MATRIX -c $ID -o pastml.ML.out --work_dir work_pastml --parameter $VAR
pastml --threads $threads -t ${tree} -d pastml.input --prediction_method MAP -m CUSTOM_RATES --rate_matrix $MATRIX -c $ID -o pastml.ML.out --work_dir work_pastml --parameter $VAR 2>&1 || true  > pastml.log

# We concatenate the aminoacids for all sites
echo ">root" > ${outprefix}_root.fasta

# Reorder by pastml numbering
head -n 1 pastml.ML.out | tr '\t' '\n' > num
# In case it's called root
cat pastml.ML.out | (grep root||:) | tr '\t' '\n' > c
# In case it's called N00000001
cat pastml.ML.out | (grep "N000000001"||:) | tr '\t' '\n' >> c
# We reorder by index + take the second columns and write it in the fasta file
paste num c | tail -n+2 | sort -k 1 -n | cut -f 2 | tr -d '\n' >> ${outprefix}_root.fasta
echo >> ${outprefix}_root.fasta

#cat pastml.ML.out | (grep root||:) | sed 's/root\t//g'| tr -d '\t' >> ${outprefix}_root.fasta
# In case it's called "N000000001"
#cat pastml.ML.out | (grep "N000000001"||:) | sed 's/N000000001\t//g'| tr -d '\t' >> ${outprefix}_root.fasta
