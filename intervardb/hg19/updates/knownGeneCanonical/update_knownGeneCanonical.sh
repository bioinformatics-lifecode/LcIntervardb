#!/bin/bash
# update_knownGeneCanonical.sh

set -e

mkdir -p ucsc_tmp && cd ucsc_tmp

wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownCanonical.txt.gz
wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz

zcat knownGene.txt.gz | awk -F'\t' '{print $1"\t"$12"\t"$8"\t"$4"\t"$5}' | sort -k1,1 > ensembl_to_ucsc.txt
zcat knownCanonical.txt.gz | awk -F'\t' '{print $5}' | sort -k1,1 > canonical_ensembl.txt

join -1 1 -2 1 -t $'\t' canonical_ensembl.txt ensembl_to_ucsc.txt | \
  awk -F'\t' '{print $2" "$3" "$4" "$5}' | sort -k1,1 -u > knownCanonical_processed.txt

echo "Transcript(ucsc/known) exons start end" > knownGeneCanonical.txt.hg19
cat knownCanonical_processed.txt >> knownGeneCanonical.txt.hg19

mv knownGeneCanonical.txt.hg19 ..
cd .. && rm -rf ucsc_tmp

echo "Done: knownGeneCanonical.txt.hg19 ($(wc -l < knownGeneCanonical.txt.hg19) lines)"
head -5 knownGeneCanonical.txt.hg19
