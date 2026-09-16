
# Getting your hands dirty: writing new Fortran code

This second part of the tutorial consists in a set of guided exercises that exemplify how to make basic additions to the code, taking as starting point the template-like, relatively self-contained submodule `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90`, which is also a good display of ALF's internal workings.

These worked-out exercises, together with ALF's modularity, boosted by its Predefined Structures, should make getting your hands dirty less daunting than it may sound.

## Downloading and using the code and tutorial

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

### Exercise: Dimensional crossover

Here we will modify the code so as to allow for different hopping matrix elements along the $x$ and $y$ directions of a square lattice.

#### Modifying the hopping

To do so we start from the submodule `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90`, which we here shorten to "`Vanilla`", found in `$ALF_DIR/Prog/Hamiltonians/`, proceeding as follows:

- Add `Ham_Ty` to the `VAR_Hubbard_Plain_Vanilla` name space in the parameter file `parameters`.
- Declare a new variable, `Ham_Ty`, in the module's specification on a new line between `!#PARAMETERS START# VAR_Hubbard_Plain_Vanilla` and the corresponding `!#PARAMETERS END#`.
- Modify the hopping matrix in the subroutine `Ham_Hop` in `Vanilla`:

```fortran
Do I = 1,Latt%N
   Ix = Latt%nnlist(I,1,0)
   Op_T(1,nf)%O(I,  Ix) = cmplx(-Ham_T,    0.d0, kind(0.D0))
   Op_T(1,nf)%O(Ix, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
   If ( L2 > 1 ) then
      Iy = Latt%nnlist(I,0,1)
!!!!!!! Modifications for Exercise 1a
      !Op_T(1,nf)%O(I,  Iy) = cmplx(-Ham_T,    0.d0, kind(0.D0))
      !Op_T(1,nf)%O(Iy, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
      Op_T(1,nf)%O(I,  Iy) = cmplx(-Ham_Ty,    0.d0, kind(0.D0))
      Op_T(1,nf)%O(Iy, I ) = cmplx(-Ham_Ty,    0.d0, kind(0.D0))
!!!!!!!
   endif
   Op_T(1,nf)%O(I,  I ) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
   Op_T(1,nf)%P(i) = i 
Enddo
```

<!-- Note: If you'd like to run the simulation using MPI, you should also add the broadcasting call for `Ham_Ty` to `Ham_Set`. -->

It is a good idea as well to get the new simulation parameter written into the file `info`, also a change in `Ham_Set`.

In the directory `Solutions/Exercise_1` we have the modified and original submodules, as well as reference data and the necessary `Start` directory (remember to copy its contents to every new `Run` directory, and to have a different `Run` directory for each simulation).

As an application of this code, we can once again consider a ladder system (e.g, a 2-leg ladder with `L1=14` and `L2=2`), for different values of `Ham_Ty`. The results you should obtain for the total spin correlation function (file `SpinT_eqJR`) are summarized in [](#fig-ladder).

```{figure} Figures/ladder.pdf
---
width: 60%
name: fig-ladder
---
Spin correlation functions along one leg for the Hubbard ladder. As $t_y$ grows the spin gap becomes large enough so as to detect the exponential decay of the spin correlation function on this small lattice size. The underlying physics of odd-even ladder systems is introduced in the article: Elbio Dagotto and T. M. Rice, "Surprises on the way from one- to two-dimensional quantum magnets: The ladder materials", Science **271**, 5249 (1996) pp. 618–623.
```

#### The SU(2) Hubbard-Stratonovich transformation

The SU(2) Hubbard-Stratonovich decomposition couples to the density and conserves spin rotational symmetry. Introduce into the module `Hamiltonian_Hubbard_smod.F90` and into the name space `VAR_Hubbard` the same changes done to the `Vanilla` module, described in the previous item, and enter `Mz=.F.` in the `parameters` file (in order to choose the $SU(N)$ Hubbard interaction) and compare results to those of the $M_z$ decomposition above — especially with regard to numerical convergence.

### Exercise: Defining a new model — The one-dimensional t-V model

In this section, we will show what modifications have to be carried out for computing the physics of the one dimensional t-V model of spinless fermions.

```{math}
\hat{H} =  -t \sum_{i} \left( \hat{c}^{\dagger}_{i}  \hat{c}^{\phantom\dagger}_{i+a}      + \hat{c}^{\dagger}_{i+a}  \hat{c}^{\phantom\dagger}_{i} \right) - \frac{V}{2} \sum_{i}  \left( \hat{c}^{\dagger}_{i}  \hat{c}^{\phantom\dagger}_{i+a}      + \hat{c}^{\dagger}_{i+a}  \hat{c}^{\phantom\dagger}_{i} \right)^2 .
```

The above form is readily included in the ALF since the interaction is written in terms of a perfect square. Expanding the square yields (up to a constant) the desired model:

```{math}
\hat{H} =  -t \sum_{i} \left( \hat{c}^{\dagger}_{i}  \hat{c}^{\phantom\dagger}_{i+a}      + \hat{c}^{\dagger}_{i+a}  \hat{c}^{\phantom\dagger}_{i} \right)  + V  \sum_{i}  \left( \hat{n}_{i} - 1/2 \right)  \left( \hat{n}_{i+a} - 1/2 \right).
```

Note that the t-V model is already implemented in ALF in the module `Hamiltonian_tV_smod.F90` found in `$ALF_DIR/Prog/Hamiltonians/`. While it can be used for checking your own results[^tV_notebook], you are not supposed to reproduce that implementation — which is more general and makes use of predefined structures — but instead to write a simpler one, based on the module `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90` as detailed below.

[^tV_notebook]: A short simulation of the t-V model can be conveniently run using its Jupyter notebook available in pyALF.

#### Define new model

In the directory `Solutions/Exercise_2` we have duplicated the ALF and commented the changes that have to be carried out to the file `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90`, which we here shorten to "`Vanilla`", found in `$ALF_DIR/Prog/Hamiltonians/`. The following are the essential steps to be carried out:

- Add the `VAR_t_V` name space in the file `parameters` and set the necessary variables — or simply rename the `VAR_Hubbard_Plain_Vanilla` name space to `VAR_t_V` and, within it, `Ham_U` to `Ham_Vint`. (Ignore the name space `VAR_tV`, which is used by the general implementation mentioned above.)
- Declare a new variable, `Ham_Vint`, in `Vanilla`'s specification.
- Add the `VAR_t_V` name space (`namelist`) declaration at the `Ham_Set` subroutine of `Vanilla`, containing the same variables the name space contains in `parameters`, and read it in.
- Still in the `Ham_set` subroutine of `Vanilla`: set `NF=1`, since we are working with spinless fermions; change the output to the info file.
- In the `Ham_V` subroutine you have to add the new interaction. For a given bond at a given time slice, we need to decouple the interaction:

```{math}
e^{\Delta \tau \frac{V}{2} \left( \hat{c}^{\dagger}_{i}  \hat{c}^{\phantom\dagger}_{i+a}      + \hat{c}^{\dagger}_{i+a}  \hat{c}^{\phantom\dagger}_{i} \right)^2} = 
\sum_{l= \pm1, \pm 2} \gamma_l e^{ \sqrt{\Delta \tau \frac{V}{2}}  \eta_l   \left( \hat{c}^{\dagger}_{i}  \hat{c}^{\phantom\dagger}_{i+a}      + \hat{c}^{\dagger}_{i+a}  \hat{c}^{\phantom\dagger}_{i} \right) } =  \sum_{l= \pm1, \pm 2} \gamma_l e^{ g  \eta_l   \left( \hat{c}^{\dagger}_{i},  \hat{c}^{\dagger}_{i+a}  \right)  O 
    \left(\hat{c}^{\phantom\dagger}_{i},  \hat{c}^{\phantom\dagger}_{i+a} \right)^{T} }.
```

Here is how this translates in the code (the new integer variable, `i2`, should be declared):

```fortran
Allocate(Op_V(Ndim,N_FL))
 do nf = 1,N_FL
    do i  =  1, Ndim
       call Op_make(Op_V(i,nf),2)
    enddo
 enddo
 Do i = 1, Ndim  ! Runs over bonds = # of lattice sites in one-dimension. 
    i2                = Latt%nnlist(i,1,0)
    Op_V(i,nf)%P(1)   = i
    Op_V(i,nf)%P(2)   = i2
    Op_V(i,nf)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.d0))
    Op_V(i,nf)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.d0)) 
    Op_V(i,nf)%g      = sqrt(cmplx(Dtau*Ham_Vint/2.d0, 0.d0, kind(0.d0)))  
    Op_V(i,nf)%alpha  = cmplx(0d0  ,0.d0, kind(0.d0))
    Op_V(i,nf)%type   = 2
    Call Op_set( Op_V(i,nf) )
 enddo
```

- Finally, you have to update the `Obser` and `ObserT` routines for the calculation of equal- and time-displaced correlations. For the `t_V` model you can essentially use the same observables as for the `Hubbard_SU(2)` model in 1D — a step which requires a number of changes with respect to the `Vanilla` base, such as:

```fortran
!!!!! Modifications for Exercise 2
 !Zpot = Zpot*ham_U                       ! Vanilla
 Zpot = Zpot*Ham_Vint                     ! t-V
 !!!!!
```

and

```fortran
!Zrho = Zrho + Grc(i,i,1) +  Grc(i,i,2)  ! Vanilla
 Zrho = Zrho + Grc(i,i,1)                 ! t-V
```

with the observables being coded in the routine `Obser` as

```fortran
Z = cmplx(dble(N_SUN), 0.d0, kind(0.D0))
 Do I1 = 1,Ndim
    I = I1
    no_I = 1
    Do J1 = 1,Ndim
       J = J1
       no_J = 1
       imj = latt%imj(I,J)
       Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
       &               Z * GRC(I1,J1,1) * ZP*ZS  ! Green
       Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
       &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS  ! SpinZ
       Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
       &               ( GRC(I1,I1,1) * GRC(J1,J1,1) * Z + &
       &                 GRC(I1,J1,1) * GR(I1,J1,1 )       ) * Z * ZP*ZS  ! Den
    enddo
    Obs_eq(3)%Obs_Latt0(no_I) =  Obs_eq(3)%Obs_Latt0(no_I) + Z * GRC(I1,I1,1) * ZP * ZS
 enddo
```

among other changes — with similar ones in the `ObserT` routine.

All necessary changes are implemented and clearly indicated in the solution provided in `Solutions/Exercise_2/Hamiltonian_Hubbard_Plain_Vanilla_smod-Exercise_2.F90`.

In the directory `Solutions/Exercise_2` we have the modified and original submodules, as well as reference data and the necessary `Start` directory (remember to copy its contents to every new `Run` directory, and to have a different `Run` directory for each simulation).

#### Phase transition

You can now run the code for various values of $V/t$. A Jordan-Wigner transformation will map the `t_V` model onto the XXZ chain:

```{math}
\hat{H} = J_{xx} \sum_{i}   \hat{S}^{x}_i \hat{S}^{x}_{i+a} +   \hat{S}^{y}_i \hat{S}^{y}_{i+a}  + J_{zz}  \sum_{i}\hat{S}^{z}_i \hat{S}^{z}_{i +a}\quad,
```

with $J_{zz} = V$ and $J_{xx} = 2t$. Hence, when $V/t = 2$ we reproduce the Heisenberg model. For $V/t > 2$ the model is in the Ising regime with long-range charge density wave order and is an insulator. In the regime $-2 < V/t < 2$ the model is metallic and corresponds to a Luttinger liquid. Finally, at $V/t < -2$ phase separation between hole rich and electron rich phases occur. [](#fig-tV) shows typical results for the density-density correlation function (file `Den_eqJR`).

```{figure} Figures/tV.pdf
---
width: 60%
name: fig-tV
---
Density-Density correlation functions of the t-V model. In the Luttinger liquid phase, $-2 < V/t < 2$, it is known that the density-density correlations decay as $\langle n(r) n(0)\rangle \propto \cos(\pi r) r^{-\left(1+K_\rho \right)}$ with $\left(1+K_\rho \right)^{-1}= \frac{1}{2} + \frac{1}{\pi} \arcsin \left( \frac{V}{2 |t|}\right)$ (A. Luther and I. Peschel, Calculation of critical exponents in two dimensions from quantum field theory in one dimension, Phys. Rev. B **12** (1975), pp. 3908–3917.) The interested reader can try to reproduce this result.
```

<!-- Challenge: How would you use the code to carry out simulations at $V/t < 0$? -->

### Exercise: Adding a new observable

This exercise illustrates the modifications that are required to implement a new observable, the correlation function of the bond-hopping in the 1-dimensional Hubbard chain. This observable is interesting in many different setups. For example, the model studied here exhibits an emergent $SO(4)$ symmetry that relates the anti-ferromagnetic order parameter and the bond dimerization (I. Affleck, PRL **55** (1985) pp. 1355; I. Affleck and F. D. M. Haldane, PRB **36** (1987) pp. 5291). Another example where this quantity is useful to investigate is the 1-dimensional Su-Schrieffer-Heeger model describing an electron-phonon system.

#### Applying Wick's theorem

Here the task is to define a new equal-time observable, the kinetic energy correlation, given by

```{math}
:label: eq-new-obs
\left\langle \hat{O}_{i,\delta} \hat{O}_{j,\delta'} \right\rangle  -  \left\langle \hat{O}_{i,\delta} \right\rangle \left\langle \hat{O}_{j,\delta'} \right\rangle  =  S_O\big(i-j,\delta,\delta'\big)
```

where $i,j$ refer to the unit cells and $\delta$ encodes the bond label. Since we are working in 1D, there is only one bond per unit cell and $\delta=ax$ such that

```{math}
\hat{O}_{i,x} = \sum_{\sigma}\left( \hat{c}^\dagger_{i,\sigma}\hat{c}^{\phantom\dagger}_{i+ax,\sigma} +\mathrm{H.c.} \right).
```

Note that the first term of Eq. {eq}`eq-new-obs` is of the generic form $\sum_{\sigma,\sigma'}\left\langle \hat{c}^\dagger_{i_1,\sigma}\hat{c}^{\phantom\dagger}_{i_2,\sigma} \hat{c}^\dagger_{j_1,\sigma'}\hat{c}^{\phantom\dagger}_{j_2,\sigma'} \right\rangle$. This expectation value can be readily decomposed into single-particle Green functions by using Wick's theorem. It can be applied for a fixed field configuration $\Phi$ since the Hamiltonian $\hat{H}(\Phi)$ is then bi-linear in the fermion operators.

```{math}
\sum_{\sigma,\sigma'}\left\langle \hat{c}^\dagger_{i_1,\sigma}\hat{c}^{\phantom\dagger}_{i_2,\sigma} \hat{c}^\dagger_{j_1,\sigma'}\hat{c}^{\phantom\dagger}_{j_2,\sigma'} \right\rangle_\Phi = 
\sum_{\sigma,\sigma'}\left(\left\langle \hat{c}^\dagger_{i_1,\sigma}\hat{c}^{\phantom\dagger}_{i_2,\sigma}\right\rangle_\Phi
\left\langle\hat{c}^\dagger_{j_1,\sigma'}\hat{c}^{\phantom\dagger}_{j_2,\sigma'} \right\rangle_\Phi
+\left\langle \hat{c}^\dagger_{i_1,\sigma}\hat{c}^{\phantom\dagger}_{j_2,\sigma'}\right\rangle_\Phi
\left\langle\hat{c}^{\phantom\dagger}_{i_2,\sigma}\hat{c}^\dagger_{j_1,\sigma'} \right\rangle_\Phi
\right)
```

The second term vanishes for $\sigma\neq\sigma'$ due to flavor symmetry (Mz-decoupling used here in the vanilla version) or due to the $SU(2)$ symmetry (density decoupling available in the generic implementation of Hubbard model). The single-particle Green functions are provided in the `Obser` routine, where all equal-time observables are measured, as

```{math}
GRC(i,j,\sigma) & = \left\langle \hat{c}^\dagger_{i,\sigma}\hat{c}^{\phantom\dagger}_{j,\sigma}\right\rangle_\Phi\\
GR(i,j,\sigma) & = \left\langle \hat{c}^{\phantom\dagger}_{i,\sigma}\hat{c}^\dagger_{j,\sigma}\right\rangle_\Phi .
```

#### Necessary code modifications

In the directory `Solutions/Exercise_3` we have the modified and original submodule `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90` found in `$ALF_DIR/Prog/Hamiltonians` (which we here shorten to "`Vanilla`"), as well as reference data and the necessary `Start` directory (remember to copy its contents to every new `Run` directory, and to have a different `Run` directory for each simulation). The following are the essential steps to be carried out:

- Introduce the new observable and allocate the memory required to store the measurements. This is done in the subroutine `Alloc_obs(Ltau)` by increasing the length of the array `Obs_eq` appropriately and adding a new case to specify the filename in which the results are stored on disc. (You might want to revisit this section later on to add the time-displaced version of the correlation function by changing `Obs_tau` in the same fashion.)
- The actual measurements are taken in the subroutine `Obser(GR,Phase,Ntau)`. While `GR` is passed to the subroutine, the first lines of code already implement the construct $\mathtt{GRC}=1-\mathtt{GR}^T$. (This section does not have to be modified, but it is useful to keep this in mind for future reference when you implement a new model from scratch.)
- The measurement of an equal-time correlation function consists of two separate parts: the connected one given by the first term in Eq. {eq}`eq-new-obs`, and the background, given by the second term.
- Implement the measurement of the connected part, stored as $\mathtt{Obs\_eq}(6)$ in this example, using the Wick decomposition sketched above.
- Keep in mind that the background $\sum_i\left\langle \hat{O}_{i,\delta} \right\rangle\neq0$ is non-vanishing and has to be measured separately (you can compare with the density correlation function), and is stored in $\mathtt{Obs\_eq}(6)\%\mathtt{Obs\_Latt0}(1)$.
- The analysis tool will then automatically combine both contributions and evaluate Eq. {eq}`eq-new-obs` using the jackknife method to estimate the mean and error of the correlation function.

The 1-D Hubbard exhibits an emergent $SO(4)$ symmetry:

```{math}
\left\langle \bar{S}(r)S(0) \right\rangle  &\sim  \frac{(-1)^r}{r}\ln^d(r)\\
\left\langle \hat{O}_{r,x} \hat{O}_{0,x} \right\rangle  -  \left\langle \hat{O}_{r,x} \right\rangle \left\langle \hat{O}_{0,x} \right\rangle  &\sim   \frac{(-1)^r}{r}\ln^\beta(r)
```

where $d=1/2$ and $\beta=-3/2$ (T. Sato, M. Hohenadler et al., Phys. Rev. B **104** (2021) L161105) and this exercise provides all the tools required to study it. Beware of the large system sizes, and therefore long run times, that are required to extract the logarithmic scaling corrections. More details are discussed in the appendix of above reference.

Finally, it is straightforward to implement the time-displaced version of this correlation function, following essentially the same steps as described above. The observable is now stored in `Obs_tau`, the measurements are taken in `ObserT`, and you can find the implementation in the solution to this exercise as well.

---

## Appendix: Excluded exercise item from "Dimensional crossover"

When running the code for the above ladder system and analyzing the spin correlation functions, `SpinZ_eqJR` and `SpinXY_eqJR` one will notice that it is hard to restore the SU(2) spin symmetry and that very long runs are required to obtain the desired equality

```{math}
\langle S^{z}_{i} S^{z}_{j}  \rangle  = \langle S^{y}_{i} S^{y}_{j}  \rangle  = \langle S^{x}_{i} S^{x}_{j}  \rangle. 
```

The structure of the output files `SpinZ_eqJR` and `SpinXY_eqJR` is described in the documentation. For the Mz Hubbard-Stratonovich transformation it is hence better to consider the improved estimator

```{math}
\langle \vec{S}_{i}  \cdot \vec{S}_{j}  \rangle
```

to compute the spin-spin correlations.

Here the aim is to include the new equal-time observable $\langle \vec{S}_{i} \cdot \vec{S}_{j} \rangle$ in the `Hubbard_Mz` code. To achieve this, you will have to carry out the following steps.

- In the subroutine `Alloc_obs` in the `Hamiltonian_example.f90` file you will have to add a new equal time observable with a call to `Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)` with `Ns = Latt%N; No = Norb; Filename ="SpinT", Nt=1, I=5`
- In the subroutine `Obser` you will have to add the Wick decomposition of this observable.

In the program `Hamiltonian_Examples.f90` to be found in the directory `Solutions/Exercise_2/Prog/` we have commented the changes that have to be carried out to add this observable. The new variable takes the name SpinT and the results you should obtain are summarized in [](#fig-ladder).



