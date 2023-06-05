#/bin/bash

# Usage:  get_content_from_notebook_latex.sh  <notebook_tex_file>

file=$1                                         # 1st argument -> "file", notebook latex filename
newfile="${file%.tex}-content.tex"              # "file" minus ".tex" extension followed by "-content.tex"

sed -e '/maketitle/,/end{document}/!d'        \
    -e '/maketitle/d;/end{document}/d'        \
    -e "s/{\([^{]*.png}\)/{${file%.tex}-\1/g" \
    -e 's/adjustimage{max size={0.9/adjustimage{max size={0.6/g' \
    $file  >  $newfile                          # Where: 
                                                # 1st: delete lines NOT in the range between expressions
                                                # 2nd: delete lines containing the expressions
                                                # 3rd: change {example.png} to {notebookname-example.png} ,
                                                #      \1 refers to the previous subexpression enclosed by
                                                #      \( \), and [^{]* matches the longest sequence of
                                                #      characters not including {
                                                # 4th: change image max size from 0.9 to 0.6
