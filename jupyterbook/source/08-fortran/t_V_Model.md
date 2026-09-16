
# Exercise: Defining a new model — The one-dimensional t-V model

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

```{figure} ../Figures/tV.pdf
---
width: 60%
name: fig-tV
---
Density-Density correlation functions of the t-V model. In the Luttinger liquid phase, $-2 < V/t < 2$, it is known that the density-density correlations decay as $\langle n(r) n(0)\rangle \propto \cos(\pi r) r^{-\left(1+K_\rho \right)}$ with $\left(1+K_\rho \right)^{-1}= \frac{1}{2} + \frac{1}{\pi} \arcsin \left( \frac{V}{2 |t|}\right)$ (A. Luther and I. Peschel, Calculation of critical exponents in two dimensions from quantum field theory in one dimension, Phys. Rev. B **12** (1975), pp. 3908–3917.) The interested reader can try to reproduce this result.
```

<!-- Challenge: How would you use the code to carry out simulations at $V/t < 0$? -->
