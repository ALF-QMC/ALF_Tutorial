#/bin/bash

file=$1                                          # 1st argument: "file"
newfile="${file%.tex}-content.tex"               # "file" minus ".tex" extension followed by "-content.tex"

sed -e '/maketitle/,/end{document}/!d' \
    -e '/maketitle/d;/end{document}/d' \
    $file  >  $newfile                           # Where: 
                                                 # 1st: delete lines NOT in the range between expressions
                                                 # 2nd: delete lines containing the expressions
