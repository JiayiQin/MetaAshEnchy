#########################################################################
# File Name: 4_build_OTU_table.sh
# Author: Jiayi Qin
# mail: qinwaii@gmail.com
# Created Time: Tue Jul  5 15:31:25 2016
#########################################################################
#!/bin/bash


FASTA="OUTPUT/lib_262/Enchy_262.fas"
SCRIPT="OTU_contingency_table.py"
STATS="${FASTA/.fas/_1f.stats}"
SWARMS="${FASTA/.fas/_1f.swarms}"
REPRESENTATIVES="${FASTA/.fas/_1f_representatives.fas}"
UCHIME="${FASTA/.fas/_1f_representatives.uchime}"
ASSIGNMENTS="${FASTA/.fas/_1f_representatives.results}"
QUALITY="OUTPUT/lib_262/Enchy.qual"
OTU_TABLE="${FASTA/.fas/_1f_OTU.table}"

python \
"${SCRIPT}" \
"${REPRESENTATIVES}" \
"${STATS}" \
"${SWARMS}" \
"${UCHIME}" \
"${QUALITY}" \
"${ASSIGNMENTS}" \
OUTPUT/lib_262/*S*.fas > "${OTU_TABLE}"
