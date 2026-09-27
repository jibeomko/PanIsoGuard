#!/usr/bin/env bash
# Rebuild every file in notes/data/ from make_toy.py. The outputs are committed, so the
# notes run without STAR or SQANTI3; this script only shows where they came from.
#
# Needs python3, samtools, STAR and SQANTI3 on PATH (I used the conda env that SQANTI3
# 6.0.1 installs: Python 3.11, samtools 1.23.1, STAR 2.7.11b).
#
#     cd notes/data && bash build_toy.sh
set -euo pipefail
cd "$(dirname "$0")"

python3 make_toy.py
samtools faidx toy.fa && samtools faidx hap1.fa && samtools faidx hap2.fa

# Long reads -> sorted, indexed BAM (the mapping axis reads FLAG/MAPQ/CIGAR only).
samtools view -b -o long_reads.bam long_reads.sam && samtools index long_reads.bam

# Short reads -> STAR with the reference annotation -> SJ.out.tab (default filters).
rm -rf star_idx star_out && mkdir -p star_idx star_out
STAR --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles toy.fa \
     --sjdbGTFfile reference.gtf --sjdbOverhang 99 --genomeSAindexNbases 4 \
     --outFileNamePrefix star_idx/ > /dev/null
STAR --genomeDir star_idx --readFilesIn short_reads.fq --outSAMtype None \
     --outFileNamePrefix star_out/ > /dev/null
cp star_out/SJ.out.tab short_reads.SJ.out.tab

# SQANTI3 QC of caller 1 against the reference annotation.
rm -rf sqanti_out
sqanti3_qc.py --isoforms caller1.gtf --refGTF reference.gtf --refFasta toy.fa \
    -o caller1 -d sqanti_out --report skip -t 1 > sqanti_out.log 2>&1
cp sqanti_out/caller1_classification.txt caller1_classification.txt

rm -rf star_idx star_out sqanti_out sqanti_out.log Log.out
echo "built: $(ls | tr '\n' ' ')"
