#########################################################################
# File Name: 2_demulti_primercut_derep_quality.sh
# Author: Jiayi Qin
# mail: qinwaii@gmail.com
# Created Time: Tue Jul  5 15:27:06 2016
# Edit:(Tue Jul 19) based on Mahe's original version
#	-control the cutadapt quality 
#		 primer anchored, output length >158, quality -e 0.15 -O 10 
#	-four library will be processed together
# Edit:(Mon Jul 25) used for Nextera library preparation result
# Updade:(Wed Sep 28) error about producing quality file was corrected.
# Upadte:(Thu Nov 10, 2016) The primer changed to H3a 
#			    Minilength was calculated according to the reference from Rüdiger, which is 262 bp
#########################################################################
#!/bin/bash

# Define variables
for f in `ls *R1_001.fastq | sed 's/_L001_R1_001.fastq//'|sed 's/JJ-//'` ; do

INPUT="OUTPUT/${f}.fastq"
PRIMER_F="AGCAGCTGGCNACVAAGGC"
PRIMER_R="ATYCACGCCAAGCGYGTCA"
MIN_LENGTH_PRIMER=262
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))  # primer match is >= 2/3 of primer length
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))

# Define binaries, temporary files and output files
CUTADAPTPRIMER="$(which cutadapt) -e 0.15 -O 10 --discard-untrimmed --match-read-wildcards --minimum-length ${MIN_LENGTH_PRIMER}"
VSEARCH=$(which vsearch)
INPUT_REVCOMP=$(mktemp)
TMP_FASTQ=$(mktemp)
TMP_FASTQ2=$(mktemp)
TMP_FASTA=$(mktemp)
OUTPUT="OUTPUT/lib_262/Enchy.output" # output couldn't define as a temperary file in this loop, if it is a permanent file, it will be kept, temperary file will be wiped out in each loop.
QUALITY_FILE="OUTPUT/lib_262/Enchy.qual"

# Reverse complement fastq file
"${VSEARCH}" --quiet \
--fastx_revcomp "${INPUT}" \
--fastqout "${INPUT_REVCOMP}"

LOG="OUTPUT/lib_262/${f}.log"
FINAL_FASTA="OUTPUT/lib_262/${f}.fas"

# Trim tags, forward & reverse primers (search normal and antisens)
cat "${INPUT}" "${INPUT_REVCOMP}" | \
     ${CUTADAPTPRIMER} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
     ${CUTADAPTPRIMER} -a "${PRIMER_R}" -O "${MIN_F}" - 2>> "${LOG}" > "${TMP_FASTQ}"

# Discard erroneous sequences and add expected error rates
"${VSEARCH}" \
--quiet \
--fastq_filter "${TMP_FASTQ}" \
--fastq_maxns 0 \
--relabel_sha1 \
--eeout \
--fastqout "${TMP_FASTQ2}" 2>> "${LOG}"

# Convert fastq to fasta (discard sequences containing Ns)
"${VSEARCH}" \
--quiet \
--fastq_filter "${TMP_FASTQ}" \
--fastq_maxns 0 \
--fastaout "${TMP_FASTA}" 2>> "${LOG}"

# Dereplicate at the study level (vsearch)
"${VSEARCH}" \
--quiet \
--derep_fulllength "${TMP_FASTA}" \
--sizeout \
--fasta_width 0 \
--relabel_sha1 \
--output "${FINAL_FASTA}" 2>> "${LOG}"

# Discard quality lines, extract sha1, expected error rates and read length
sed 'n;n;N;d' "${TMP_FASTQ2}" | \
awk 'BEGIN {FS = "[;=]"}
{if (/^@/) {printf "%s\t%s\t", $1, $3} else {print length($1)}}' | \
tr -d "@" >> "${OUTPUT}"

done 

# Produce the final quality file
sort -k3,3n -k1,1d -k2,2n "${OUTPUT}" | \
uniq --check-chars=40 > "${QUALITY_FILE}"

# Clean
rm -f "${INPUT_REVCOMP}" "${TMP_FASTQ}" "${TMP_FASTA}" "${TMP_FASTQ2}" #"${OUTPUT}"
