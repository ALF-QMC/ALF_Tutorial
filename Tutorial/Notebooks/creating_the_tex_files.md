
1) Download latex version of the notebook from Jupyter server (in the browser):
```
example.tex
```

2) Apply `get_content_from_notebook_latex.sh` script:
```
get_content_from_notebook_latex.sh example.tex
```
which outputs a file called `example-content.tex` with only the lines of `example.tex` found between `\maketitle` and `\end{document}`.


---

Content of script `get_content_from_notebook_latex.sh`:

```bash
#/bin/bash

file=$1                                          # 1st argument: "file"
newfile="${file%.tex}-content.tex"               # "file" minus ".tex" extension followed by "-content.tex"

sed -e '/maketitle/,/end{document}/!d' \
    -e '/maketitle/d;/end{document}/d' \
    $file  >  $newfile                           # Where: 
                                                 # 1st: delete lines NOT in the range between expressions
                                                 # 2nd: delete lines containing the expressions

```

---

* The file `notebooks_preamble.tex` is the preamble of the downloaded example.tex, minus the lines:
```
\documentclass[11pt]{article}
\usepackage{mathpazo}
```
