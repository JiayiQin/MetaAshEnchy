#########################################################################
# File Name: 1_merge.sh
# Author: Jiayi Qin
# mail: qinwaii@gmail.com
# Created Time: Tue Jul  5 15:24:30 2016
#########################################################################
#!/bin/bash
for f in `ls *R1_001.fastq | sed 's/_L001_R1_001.fastq//'|sed 's/JJ-//'` ; do

VSEARCH=$(which vsearch)
THREADS=4
ENCODING=33
FORWARD="JJ-${f}_L001_R1_001.fastq"  # bz2, gz and uncompressed fastq files are allowed
REVERSE="JJ-${f}_L001_R2_001.fastq"
FORWARD_TRIM=$(mktemp)
REVERSE_TRIM=$(mktemp)
SINGLE_TRIM=$(mktemp)
OUTPUT="OUTPUT/$f.fastq"

#module load sickle/1.33

#Remove bad 3' stretch
sickle pe -x -f ${FORWARD} -r ${REVERSE} \
		-t sanger -o ${FORWARD_TRIM} 2>${OUTPUT/.fastq/.log} -p ${REVERSE_TRIM}\
		-s ${SINGLE_TRIM} 

# Merge read pairs
"${VSEARCH}" \
	--threads ${THREADS} \
	--fastq_mergepairs ${FORWARD_TRIM} \
	--reverse ${REVERSE_TRIM} \
	--fastq_ascii ${ENCODING} \
	--fastqout ${OUTPUT} \
	--fastq_allowmergestagger \
	--quiet 2>> ${OUTPUT/.fastq/.log}

rm ${FORWARD_TRIM} ${REVERSE_TRIM} 
rm ${SINGLE_TRIM}
done
