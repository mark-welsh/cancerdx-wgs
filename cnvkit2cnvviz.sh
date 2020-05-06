#!/bin/bash

CNR_FILE="$1"
B_ALLELE_VCF="$2"
NUM_CORES="$3"

cavatica_name="$(echo ${CNR_FILE} | awk -F'/' '{print $NF}' | cut -d'.' -f1 | cut -d'-' -f1)"
tumor_name="$(echo ${CNR_FILE} | awk -F'/' '{print $NF}' | cut -d'.' -f1 | cut -d'-' -f2 | tr '_' '-')"
normal_name="$(echo ${CNR_FILE} | awk -F'/' '{print $NF}' | cut -d'.' -f1 | cut -d'-' -f3 | tr '_' '-')"

basename="$(echo ${cavatica_name}"_"${tumor_name}"_WGSRD_"${cavatica_name}"_"${normal_name}"_WGSRD")"

awk '{print $1"\t"$2"\t"$3"\t"$6"\t"".""\t"$4"\t""NA""\t""NA"}' "${CNR_FILE}" > ${basename}.txt
tail -n +2 ${basename}.txt > ${basename}.bed

bedtools intersect -header -a ${B_ALLELE_VCF} -b "${basename}.bed" | \
    bedtools intersect -header -a stdin -b /usr/local/bin/resources/b_allele.expanded.bed | \
    bcftools norm -m -both | \
    bcftools query -f '%CHROM\t%POS\t[%AF]\n' | sort -u >> ${basename}.middle.tsv
 
sed -E -i "" 's/^chr(.*)/\1/' ${basename}.bed
sed -E -i "" 's/^chr(.*)/\1/' ${basename}.middle.tsv

python /usr/local/bin/make_panels.py \
    --num_cores ${NUM_CORES} \
    --top_panel ${basename}.bed \
    --mid_panel ${basename}.middle.tsv \
    --sample_id ${basename}

sed -E -i '1s/.*/#top/' ${basename}.panels.txt
sed -E -i '1s/.*/#middle/' ${basename}.middle.panels.txt
cat ${basename}.middle.panels.txt >> ${basename}.panels.txt

rm -f ${basename}.txt ${basename}.bed ${basename}.top.txt ${basename}.middle.tsv ${basename}.middle.panels.txt
