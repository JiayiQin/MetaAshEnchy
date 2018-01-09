# filter the blast output
# Wed Dec 13 19:10:23 CET 2017

INPUT="OUTPUT/lib_262/Enchy_262_1f_OTU.filtered.local.10output.262length.blast"
OUTPUT="${INPUT/.blast/.filtered.blast}"

awk '$3 > 97 && $4 >= 262' "${INPUT}" > "${OUTPUT}"

