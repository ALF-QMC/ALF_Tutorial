
Create usable Latex notebook files
---

Three steps:

1. Download the latex version of the notebook from the Jupyter server [ in the browser: File > Download as > LaTeX (.tex) ], say `example01.zip`, which contains:
```
example01.tex
output_12_2.png
```

2. Rename any figures by preceding their original name with `example01-`, i.e., `$(basename <filename> .png)-`:
```
mv output_12_2.png example01-output_12_2.png 
```

3. Apply `get_content_from_notebook_latex.sh` script:
```
get_content_from_notebook_latex.sh example01.tex
```
which outputs a file called `example01-content.tex` with only the lines of `example01.tex` found between `\maketitle` and `\end{document}`. It also changes `{example.png}` to `{example01-example.png}`, since different downloaded latex notebooks may have figures with the same name.


---

Remarks
---

- The file **`notebooks_preamble.tex`** is the preamble of the downloaded `minimal_ALF_run.tex`, minus the lines:
```
\documentclass[11pt]{article}
\usepackage{mathpazo}
\usepackage{titling}
\DeclareCaptionFormat{nocaption}{}
\captionsetup{format=nocaption,aboveskip=0pt,belowskip=0pt}
```

- TO DO: instead of downloading notebooks' html and latex versions, use nbconvert command line

- TO DO: include in the pipeline the creation of notebooks' html and latex versions


---

Content of script `get_content_from_notebook_latex.sh`:
---
```bash

#/bin/bash

# Usage:  get_content_from_notebook_latex.sh  <notebook_tex_file>

file=$1                                         # 1st argument -> "file", notebook latex filename
newfile="${file%.tex}-content.tex"              # "file" minus ".tex" extension followed by "-content.tex"

sed -e '/maketitle/,/end{document}/!d'        \
    -e '/maketitle/d;/end{document}/d'        \
    -e "s/{\([^{]*.png}\)/{${file%.tex}-\1/g" \
    $file  >  $newfile                          # Where: 
                                                # 1st: delete lines NOT in the range between expressions
                                                # 2nd: delete lines containing the expressions
                                                # 3rd: change {example.png} to {notebookname-example.png} ,
                                                #      \1 refers to the previous subexpression enclosed by
                                                #      \( \), and [^{]* matches the longest sequence of
                                                #      characters not including {

```
