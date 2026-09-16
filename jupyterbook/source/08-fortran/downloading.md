
# Downloading and using the code and tutorial

One can use the ALF package downloaded automatically by the Python script in the first part of this tutorial, or manually, by typing

```bash
git clone https://git.physik.uni-wuerzburg.de/ALF/ALF.git
```

in a shell. Similarly, to download the tutorial, including solutions, enter:

```bash
git clone https://git.physik.uni-wuerzburg.de/ALF/ALF_Tutorial.git
```

This document is found in the directory `Tutorial-ALF-{tutALFver}`.

The necessary environment variables and the directives for compiling the code are set by the script `configure.sh`:

```bash
source configure.sh GNU noMPI
```

followed by the command `make`. Details and further options are described in the package's documentation found in its repository.

A workflow you can adopt for solving the exercises — or indeed using ALF in general — is the following:

1. Compile the modified Hamiltonian module, for instance:  
   `make all`
2. Create a data directory with the content of `Start`:  
   `cp -r ./Scripts_and_Parameters_files/Start ./Run  &&  cd ./Run/`
3. Run its executable, e.g., serially:  
   `$ALF_DIR/Prog/ALF.out`
4. Perform default analyses[^analysis_sh]:  
   `$ALF_DIR/Analysis/ana.out *`

[^analysis_sh]: The `analysis.sh` bash script from earlier versions of ALF, run without arguments, is still available in the `Start` directory.

The structure of the data files and details on the analysis output can be found in ALF's documentation.

<!-- Two common pitfalls are:
- forgetting to run `configure.sh`, and
- not using different `Run` data directories for independent runs.
-->
