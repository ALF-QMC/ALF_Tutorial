# Installation

This section shows briefly how to set your machine up to work through the rest of this tutorial for more detail, refer to the documentations of [ALF](https://alf.physik.uni-wuerzburg.de/doc.pdf) and [pyALF](https://alf.physik.uni-wuerzburg.de/pyalf-doc/).

## ALF prerequisites
To install the ALF dependencies.

:::::{dropdown} **Linux**
Execute one of these installation command on your shell.

::::{tab-set}
:::{tab-item} Debian/Ubuntu/Linux Mint
:sync: deb
```sh
sudo apt install make gfortran libopenblas-dev \
           python3 libopenmpi-dev g++ curl libghc-zlib-dev \
           git ca-certificates cmake bash
```
:::
:::{tab-item} Red Hat/Fedora/CentOS/Rocky Linux/AlmaLinux
:sync: rh
```sh
sudo dnf install gcc-gfortran make lapack-devel \
           python3 openmpi-devel gcc-c++ curl zlib-ng-compat-devel \
           git cmake 
```
:::
:::{tab-item} OpenSuSE/SLES
:sync: suse
```sh
sudo zypper install make gcc-fortran gcc-c++ lapack-devel \
             python3 mpich-devel curl zlib-devel cmake git
```
:::
:::{tab-item} Arch Linux
:sync: arch
```sh
sudo pacman -S make gcc-fortran lapack python openmpi curl zlib git
```
:::
::::
:::::

:::::{dropdown} **Other Unixes**
GNU Fortran  and the LAPACK implementation from netlib.org should be available for your system. Consult the documentation of your system on how to install the relevant packages. The package names from linux should give good starting points for your search.
:::::

:::::{dropdown} **MacOS**
For MacOS, you will need to have Xcode as well as the Apple developer tools installed. The rest of the necessary packages can e.g. be installed via [Homebrew](https://brew.sh/):
```sh
brew install gcc mpich cmake
```
:::::

:::::{dropdown} **Windows**
For Windows, we recommend using the Windows Subsystem for Linux, which can be installed following [this guide](https://learn.microsoft.com/en-us/windows/wsl/install) and then proceed according to the instructions for Linux.
:::::

## pyALF installation

```{warning}
In previous versions of pyALF, the installation instructions asked the users to set the environment variable `PYTHONPATH`.
This conflicts with the newer pip package, therefore you should remove definitions of the `PYTHONPATH` environment variable related to pyALF.
```

pyALF can be installed via the Python package installer [pip](https://pip.pypa.io/en/stable/). If you do not yet have a Python environment with working pip, we recommend [miniforge](https://conda-forge.org/download/).
Miniforge is recommended due to its convenience, but any other solution, e.g. [uv](https://docs.astral.sh/uv/) or the system package manager with Python-venv will also work.

```sh
pip install pyALF
```

## Setting ALF directory through environment variable

Since pyALF is set up to automatically clone ALF with git, it is not strictly necessary to download ALF manually, but pyALF will download ALF every time it does not find it. Therefore it is recommended to clone ALF once manually from [here](https://github.com/ALF-QMC/ALF) and setting its location in the environment variable `ALF_DIR`. This way, pyALF will use the same ALF source code directory every time.

ALF can be cloned with the Unix shell command

```bash
git clone https://github.com/ALF-QMC/ALF.git
```

This will create a folder called `ALF` in the current working directory of the terminal and download the repository there.

The environment variable can then be set with:

::::{tab-set}

:::{tab-item} Bash
:sync: bash

```bash
export ALF_DIR="/path/to/ALF"
```

where `/path/to/ALF` is the location of the ALF code, for example `/home/jonas/Programs/ALF`.

Add the line to `~/.bashrc` (Linux, interactive shells) or to
`~/.bash_profile` / `~/.profile` (macOS, or login shells), then reload with:

```bash
source ~/.bashrc
```
:::

:::{tab-item} Zsh
:sync: zsh

```zsh
export ALF_DIR="/path/to/ALF"
```

where `/path/to/ALF` is the location of the ALF code, for example `/home/jonas/Programs/ALF`.

Add the line to `~/.zshrc` (default shell on macOS), then reload with:

```zsh
source ~/.zshrc
```
:::

:::{tab-item} Fish
:sync: fish

```fish
set -gx ALF_DIR "/path/to/ALF"
```

where `/path/to/ALF` is the location of the ALF code, for example `/home/jonas/Programs/ALF`.

Add the line to `~/.config/fish/config.fish`, or set it permanently in one go
(no file editing needed):

```fish
set -Ux ALF_DIR "/path/to/ALF"
```
:::

:::{tab-item} Tcsh / Csh
:sync: tcsh

```tcsh
setenv ALF_DIR "/path/to/ALF"
```

where `/path/to/ALF` is the location of the ALF code, for example `/home/jonas/Programs/ALF`.

Add the line to `~/.tcshrc` (or `~/.cshrc`), then reload with:

```tcsh
source ~/.tcshrc
```
:::

:::{tab-item} PowerShell
:sync: powershell

```powershell
$env:ALF_DIR = "C:\path\to\ALF"
```

where `C:\path\to\ALF` is the location of the ALF code, for example `C:\Users\jonas\Programs\ALF`.

For the current user only, this lasts for the session. To persist it, add the
line to your profile (`$PROFILE`, e.g.
`Documents\PowerShell\Microsoft.PowerShell_profile.ps1`) or register it
permanently:

```powershell
[Environment]::SetEnvironmentVariable("ALF_DIR", "C:\path\to\ALF", "User")
```
:::

:::{tab-item} Windows cmd.exe
:sync: cmd

```bat
set ALF_DIR=C:\path\to\ALF
```

where `C:\path\to\ALF` is the location of the ALF code, for example `C:\Users\jonas\Programs\ALF`.

This only affects the current session. To store it permanently for your user
(takes effect in new terminals):

```bat
setx ALF_DIR "C:\path\to\ALF"
```
:::

::::

```{tip}
Check that it worked with `echo $ALF_DIR` (POSIX shells / fish),
`echo $env:ALF_DIR` (PowerShell) or `echo %ALF_DIR%` (cmd).
```

## Using Jupyter Notebooks

The tutorial heavily utilizes [Jupyter Notebooks](https://jupyter.org/) that you will have to download and execute yourself. For this, we recommend to install either [JupyterLab](https://jupyter.org/install) or [Visual Studio Code](https://code.visualstudio.com/).

## Test setup

To check if most things have been set up correctly, the script `minimal_ALF_run` can be used. It utilizes pyALF to run a short ALF simulation in the folder `ALF_data`. It ships with pyALF, one should therefore be able to run it by executing

```sh
minimal_ALF_run
```

in the Unix shell. If it does clone the ALF repository, `ALF_DIR` has not been set up correctly. Note that on the first compilation, ALF downloads and compiles HDF5, which can take several minutes.

To test the execution of Jupyter Notebooks on can download the "minimal ALF run" example from the next page through the download button on the top right (![Download button](download.svg)).
