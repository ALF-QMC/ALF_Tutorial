# For odd lines (1st, 3rd, etc.) print element in 1st column,
# for even lines (2nd, 4th, etc.) print elements in 3rd and 4th columns

# Usage example:  reformat_data_files.sh Run_V2.5/Den_eqJR

#!/bin/bash

file=$1
awk '{ printf "%s", ( NR % 2 == 1 ? $1 "\t" : $3 "\t" $4 "\n" ) }' $file > ${file%.dat}-reformat.dat
