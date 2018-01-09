
TABLE="RECHECK/lib_262/Enchy_262_1f_OTU.table"
FILTERED="${TABLE/.table/.filtered.table}"
MY_FAVORITE_TAXA="NA"
FASTA_FILE="${TABLE/.table/.filtered.fas}"

head -n 1 "${TABLE}" > "${FILTERED}"
grep "${MY_FAVORITE_TAXA}" "${TABLE}" | \
awk '$7 == "N" && $9 <= 0.0002 && ($2 >= 3 || $8 >= 2)' >> "${FILTERED}"

awk -F "\t" 'NR>1{print ">"$4"\n"$10}' "${FILTERED}">"${FASTA_FILE}"
