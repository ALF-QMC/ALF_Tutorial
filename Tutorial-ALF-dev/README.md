# ALF Tutorial
This is a tutorial for using [ALF](https://git.physik.uni-wuerzburg.de/ALF/ALF/) intended to get you started from zero up to your own first projects.

## Installation.
ALF is pretty self-contained, you only need a LAPACK and a BLAS implementation as external libraries. For compiling the source code you need make and a Fortran 2003 compatible Compiler.

pyALF, the Python interface for ALF, also demands only Python, Jupiter and a few basic Python packages: SciPy, NumPy and matplotlib.

### Dependencies

#### Python

Python and its packages can be easily installed on a variety of platforms using the Anaconda distribution -- check its [installation instructions](https://docs.anaconda.com/anaconda/install/) for your system. Then, from Anaconda, all that is needed is to issue the command
```bash
conda install -c anaconda  ipython jupyterlab scipy numpy matplotlib
```
Anaconda is recommended due to its convenience, but the system's package management (e.g., apt-get) or Python's own package management, pip3, can be used instead if preferred.

#### ALF packages

* **Linux**   
  To install the relevant ALF packages.
  - **Debian/Ubuntu/Linux Mint**:  `sudo apt-get install gfortran liblapack-dev make`
  - **Red Hat/Fedora/CentOS**:  `sudo dnf install gcc-gfortran make liblapack-devel`
  - **OpenSuSE/SLES**:  `sudo zypper install gcc-gfortran make lapack-devel`
  - **Arch Linux**:  `pacman -S make gcc-fortran lapack`

* **Other Unixes**   
  gfortran and the lapack implementation from netlib.org should be available for your system. Consult the documentation of your system on how to install the relevant packages. The package names from linux should give good starting points for your search.

* **MacOS**   
  gfortran for MacOS can be found at https://gcc.gnu.org/wiki/GFortranBinaries#MacOS. Detailed information on how to install the package  can be found at: https://gcc.gnu.org/wiki/GFortranBinariesMacOS. You will need to have Xcode as well as the  Apple developer tools installed. 

* **Windows**   
  For Windows, we recommend using the Windows Subsystem for Linux, which can be installed following this guide https://learn.microsoft.com/en-us/windows/wsl/install and then proceed accoring to the instructions for Linux.

### Building
After you have obtained the source code of [ALF](https://git.physik.uni-wuerzburg.de/ALF/ALF/) and [pyALF](https://git.physik.uni-wuerzburg.de/ALF/pyALF/), and have set up the build environment, you need to build the source files. Executing `make` in the ALF's root directory does the job, and you can check the tutorial for more details.

## Editors
For the coding parts of the exercises we recommend to use a text editor. Linux usually have one installed like kwrite, kate or gedit(emacs and VI work great, too) and for the Windows users we recommend Notepad++. Jupyter notebooks of course can be edited in a browser.
