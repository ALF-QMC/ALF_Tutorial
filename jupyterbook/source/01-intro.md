# ALF Tutorial (work in progress)

<!-- ## Todo
 1. We need  a clear reference  to the alf-doc.
 2. Where  so we include t
:::{attention} Work in progress
This is an unfinished draft of a new version of the tutorial based on [Jupyter Book](https://jupyterbook.org).
Up until [Sec. 6](08-fortran.md), the documents revolves around [Jupyter Notebooks](https://jupyter.org/) that can be downloaded to run locally. They are interactive documents that combine code, narration and display of data in one file. -->

This is a tutorial for using [ALF](https://github.com/ALF-QMC/ALF/) and [pyALF](https://github.com/ALF-QMC/pyALF/) intended to get you started from zero up to your own first projects.

## Todo   Deadline Friday 18th

- [ ] We need  a clear reference  to the alf-doc.
- [ ] Where  so we include t
- [ ] [Introduction](01-intro.md)
- [ ] [Installation](02-installation.md)  **Jonas** 
- [ ] [A minimal ALF run](03-minimal_ALF_run.ipynb) **Jonas** 
- [ ] [Reproducing some nice physics](04-nice-physics.md)  
  - [ ] [Bose-Einstein Condensate](04-nice-physics/bec.ipynb) **Luis**
  - [ ] Other ideas?  **all**
- [ ] [Monte-Carlo pitfalls and comparison to exact results](05-caveats.md)
  - [x] [Trotter systematic error - Hubbard on the square lattice](05-caveats/trotter_error.ipynb) 
  - [x] [Testing against ED - Hubbard on a ring](05-caveats/testing_against_ED.ipynb)
  - [ ] [Warmup time](05-caveats/warmup.ipynb) **Helke**+**Jonas** Use PyALF internal functions.  Choose  a observable such as spinT  at AFM wave vector.  This longer autocrrelation times as the Energy.
  - [ ] [Autocorrelation time](05-caveats/autocorrelation.ipynb)**Helke**+**Jonas** Use PyALF internal functions.  Choose  a observable such as spinT  at AFM wave vector.  This longer autocrrelation times as the Energy. 
  - [x] [Imaginary time discretization error](05-caveats/delta-tau.ipynb)  Maybe a bit redundant
- [ ] [Projector versus finite temperature](06-projector.md) **Jokob**
- [ ] [My very first sign problem](07-sign.md) **Adrien**
- [ ] [Getting your hands dirty: writing new Fortran code](08-fortran.md) **Moritz** + **Jakob**   First step is to split each exercise into a separate md file. Here we do not need  Jupyter notebooks since we have to edit code.  So this part can be directly be taken from the old .pdf.  As a first step nothing new has to be done. 

