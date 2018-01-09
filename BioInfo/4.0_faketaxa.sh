#########################################################################
# File Name: faketaxa.sh
# Author: Jiayi Qin
# mail: qinwaii@gmail.com
# Created Time: Tue Jul  5 15:30:40 2016
#########################################################################
#!/bin/bash

grep "^>" OUTPUT/lib_262/Enchy_262_1f_representatives.fas | \
sed 's/^>//
s/;size=/\t/
s/;$/\t100.0\tNA\tNA/' > OUTPUT/lib_262/Enchy_262_1f_representatives.results
