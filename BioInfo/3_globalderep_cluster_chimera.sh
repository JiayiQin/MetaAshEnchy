#########################################################################
# File Name: 3_globalderep_cluster_chimera.sh
# Author: Jiayi Qin
# mail: qinwaii@gmail.com
# Created Time: Tue Jul  5 15:29:56 2016
#########################################################################
#!/bin/bash

VSEARCH=$(which vsearch)
#SWARM=$(~/tool/swarm) please attention, this usage is not valid in Risø server
TMP_FASTA=$(mktemp --tmpdir=".")
FINAL_FASTA="OUTPUT/lib_262/Enchy_262.fas"

# Pool sequences

cat OUTPUT/lib_262/*S*.fas > "${TMP_FASTA}"

# Dereplicate (vsearch)
"${VSEARCH}" --derep_fulllength "${TMP_FASTA}" \
--sizein \
--sizeout \
--fasta_width 0 \
--output "${FINAL_FASTA}" > /dev/null

rm -f "${TMP_FASTA}"

# Clustering
THREADS=16
TMP_REPRESENTATIVES=$(mktemp --tmpdir=".")
swarm  \
-d 1 -f -t ${THREADS} -z \
-i ${FINAL_FASTA/.fas/_1f.struct} \
-s ${FINAL_FASTA/.fas/_1f.stats} \
-w ${TMP_REPRESENTATIVES} \
-o ${FINAL_FASTA/.fas/_1f.swarms} < ${FINAL_FASTA}

# Sort representatives
"${VSEARCH}" --fasta_width 0 \
--sortbysize ${TMP_REPRESENTATIVES} \
--output ${FINAL_FASTA/.fas/_1f_representatives.fas}
rm ${TMP_REPRESENTATIVES}

# Chimera checking
REPRESENTATIVES=${FINAL_FASTA/.fas/_1f_representatives.fas}
UCHIME=${REPRESENTATIVES/.fas/.uchime}
"${VSEARCH}" --uchime_denovo "${REPRESENTATIVES}" \
--uchimeout "${UCHIME}"
