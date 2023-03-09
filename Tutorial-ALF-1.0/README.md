# ALF Tutorial
This is a tutorial for using ALF intended to get you started with your first 
projects.

## Installation.
ALF is pretty self-contained, you only need a LAPACK and a BLAS implementation
as external libraries. For compiling the source code you need make and a
Fortran 2003 compatible Compiler.

### Linux
In the following we give hints on how to install relevant packages.

#### Debian/Ubuntu/Linux Mint
- sudo apt-get install gfortran liblapack-dev make

#### Red Hat/Fedora/CentOS
- sudo dnf install gcc-gfortran make liblapack-devel

#### OpenSuSE / SLES
- sudo zypper install gcc-gfortran make lapack-devel

#### Arch Linux
- pacman -S make gcc-fortran lapack


### Other Unixes
gfortran and the lapack implementation from netlib.org should be available for
your system. Consult the documentation of your system on how to install the
relevant packages. The package names from linux should give good starting points
for your search.

### MacOS
gfortran for MacOS can be found at https://gcc.gnu.org/wiki/GFortranBinaries#MacOS.
Detailed information on how to install the package  can be found at: https://gcc.gnu.org/wiki/GFortranBinariesMacOS.
You will need to have Xcode as well as the  Apple developer tools installed. 

### Windows
For Windows, we recommend using the Windows Subsystem for Linux, which can be installed following this guide https://learn.microsoft.com/en-us/windows/wsl/install and then proceed accoring to the instructions for Linux.

## Building
After you have obtained the source code and have set up the build environment
you need to build the source files. Executing make in the root directory of ALF
does the Job.

## Editors
For the coding parts of the exercises we recommend to use a text editor. 
Linux usually have one installed like kwrite, kate or gedit(emacs and VI work
great, too). For the windows users we recommend Notepad++.
