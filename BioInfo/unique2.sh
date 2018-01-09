INPUT="/Users/JiayiQin/Dropbox/2nd_Miseq/RECHECK/Enchy_262_1f_OTU.filtered.local.10output.262length.filtered.blast"
OUTPUT="${INPUT/.blast/.top.blast}"

head -n 1 "${INPUT}" > "${OUTPUT}"
grep -v "amplicon" "${INPUT}"|awk -F "\t" '!a[$1]++' >>"${OUTPUT}" #The expression select the top row for each amplicon
