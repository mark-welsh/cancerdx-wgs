#!/bin/bash

CNR_FILE="$1"
VARDICT_VCF="$2"
NUM_CORES="$3"

cavatica_name="$(echo ${CNR_FILE} | awk -F'/' '{print $NF}' | cut -d'.' -f1 | cut -d'-' -f1)"
tumor_name="$(echo ${CNR_FILE} | awk -F'/' '{print $NF}' | cut -d'.' -f1 | cut -d'-' -f2 | tr '_' '-')"
normal_name="$(echo ${CNR_FILE} | awk -F'/' '{print $NF}' | cut -d'.' -f1 | cut -d'-' -f3 | tr '_' '-')"

basename="$(echo ${cavatica_name}"_"${tumor_name}"_WGSRD_"${cavatica_name}"_"${normal_name}"_WGSRD")"

awk '{print $1"\t"$2"\t"$3"\t"$6"\t"".""\t"$4"\t""NA""\t""NA"}' "${CNR_FILE}" > ${basename}.txt
tail -n +2 ${basename}.txt > ${basename}.bed

echo "#middle" >> ${basename}.middle.txt
bedtools intersect -b "${basename}.bed" -a ${VARDICT_VCF} | \
    bcftools norm -m -both "${VARDICT_VCF}" | \
    bcftools query -f '%CHROM\t%POS\t[%AF\t]\n' | \
    awk '{print $1"\t"$2"\t"$3"\n"$1"\t"$2"\t"$4}' >> ${basename}.middle.txt
   
sed -E -i 's/^chr([0-920-22XY]\s+.*)/\1/' ${basename}.bed
sed -E -i 's/^chr([0-920-22XY]\s+.*)/\1/' ${basename}.middle.txt

python /usr/local/bin/make_panels.py \
    --num_cores ${NUM_CORES} \
    --top_panel ${basename}.bed \
    --mid_panel ${basename}.middle.txt \
    --sample_id ${basename}

sed -E -i '1i\#top' ${basename}.panels.txt
sed -E -i '1i\#middle' ${basename}.middle.panels.txt
cat ${basename}.middle.panels.txt >> ${basename}.panels.txt

rm -f ${basename}.txt ${basename}.bed ${basename}.top.txt ${basename}.middle.txt ${basename}.middle.panels.txt
