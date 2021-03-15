# ALF Tutorial
This is a tutorial for using [ALF 2.0](https://git.physik.uni-wuerzburg.de/ALF/ALF/-/tree/ALF-2.0) intended to get you started from zero up to your own first 
projects.

## Installation.
ALF is pretty self-contained, you only need a LAPACK and a BLAS implementation
as external libraries. For compiling the source code you need make and a
Fortran 2003 compatible Compiler.

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
  - **Red Hat/Fedora/CentOS**:  `sudo yum install gcc-gfortran make liblapack-devel`
  - **OpenSuSE/SLES**:  `sudo zypper install gcc-gfortran make lapack-devel`
  - **Arch Linux**:  `pacman -S make gcc-fortran lapack`

* **Other Unixes**   
  gfortran and the lapack implementation from netlib.org should be available for your system. Consult the documentation of your system on how to install the relevant packages. The package names from linux should give good starting points for your search.

* **MacOS**   
  gfortran for MacOS can be found at https://gcc.gnu.org/wiki/GFortranBinaries#MacOS. Detailed information on how to install the package  can be found at: https://gcc.gnu.org/wiki/GFortranBinariesMacOS. You will need to have Xcode as well as the  Apple developer tools installed. 

* **Windows**   
  The easiest way to compile Fortran code in Windows is trough Cygwin, which provides a Unix-like environment for Windows. The installer also works as a package manager, providing an extensive collection of software from the Unix ecosystem. For convenience, we provide a zip archive containing the Cywin installer and a local repository with all the additional software needed for 
ALF. 

  Steps for installing Cygwin:
  - Download zip from https://www.dropbox.com/s/ap8vl85gn9nfbo7/cygwin_ALF.zip?dl=0 and unzip
  - Execute "setup-x86_64.exe". If administrator rights are missing, execute it from the command line as "setup-x86.exe --no-admin".
  - In the setup choose "Install from local directory" instead of "Install from Internet".
  - Choose root directory, where cygwin will be installed. You should memorize this diretory.
  - Choose the diretory "cygwin\_ALF" (The one which also contains "setup-x86_64.exe") as local package Directory.
  - At the "Select Packages" screen, in "Categories" view, at the line marked "All", click on the word "default" so that it changes to "install".
  - Finish installation
  - Optional: To add, remove or update installed packages, rerun the installer setup and chose "Install from Internet".
  - You can now use the installed Cygwin packages by starting the Cygwin terminal. It is a UNIX terminal which, by default, starts in the home diretory "/home/<username>" of the UNIX-system emulated by Cygwin, where "/" is the root directory of Cygwin. For example, if you have installed Cygwin in "C:\cygwin64\", then the home Directory of Cygwin ca be found at "C:\cygwin64\home\<username>", in the Windows system.

### Building
After you have obtained the source code of [ALF 2.0](https://git.physik.uni-wuerzburg.de/ALF/ALF/-/tree/ALF-2.0) and the correspondent [pyALF](https://git.physik.uni-wuerzburg.de/ALF/pyALF/-/tree/ALF-2.0), and have set up the build environment, you need to build the source files. Executing `make` in the ALF's root directory does the job, and you can check the tutorial for more details.

## Editors
For the coding parts of the exercises we recommend to use a text editor. Linux usually have one installed like kwrite, kate or gedit(emacs and VI work great, too) and for the Windows users we recommend Notepad++. Jupyter notebooks of course can be edited in a browser.
