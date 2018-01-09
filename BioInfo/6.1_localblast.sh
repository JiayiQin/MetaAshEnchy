
INPUT="OUTPUT/lib_262/Enchy_262_1f_OTU.filtered.fas"

REF="/DATA_2/jiayi/Enchy_ECT/H3.fasta"

BLAST="${INPUT/.fas/.local.10output.262length.blast}"


#Blast the OTU table
blastn  -db ${REF} \
	-query ${INPUT} \
	-searchsp 262 \
	-outfmt "6 std stitle salltitles" -out ${BLAST} \
	-max_target_seqs 10 -num_threads=18 -evalue 0.001
