#########################################################################
# File Name: 0_check.sh
# Author: Jiayi Qin
# mail: qinwaii@gmail.com
# Created Time: Tue Jul  5 15:23:00 2016
#########################################################################
#!/bin/bash
VSEARCH=$(which vsearch)
FORWARD="JJ-En100_S93_L001_R1_001.fastq"  # --fastq_chars requires uncompressed fastq files

OUTPUT="OUTPUT/EnchyCODE.fastq"

# Check quality encoding (33 or 64?)
"${VSEARCH}" \
	    --fastq_chars ${FORWARD} 2> ${OUTPUT/.fastq/.log}



