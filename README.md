# ALF Tutorial
This is a tutorial for using ALF intended to get you started with your first 
projects.

## Installation.
ALF is pretty self-contained, you only need a LAPACK and a BLAS implementation
as external libraries. For compiling the source code you need make and a
Fortran03 compatible Compiler.

### Linux
In the following we give hints on how to install relevant packages.

#### Debian/Ubuntu/Linux Mint
- sudo apt-get install gfortran liblapack-dev make

#### Red Hat/Fedora/CentOS
- sudo yum install gcc-gfortran make liblapack-devel

#### OpenSuSE / SLES
- sudo zypper install gcc-gfortran make lapack-devel

ToDo:
Linux Mint: Ubuntu based
Arh Linux. Has a couple of derivatives. We will assume they are equal...

### Other Unices
gfortran and the lapack implementation from netlib.org should be available for
your system. Consult the documentation of your system on how to install the
relevant packages. The package names from linux should give good starting points
for your search.

### MacOS

### Windows