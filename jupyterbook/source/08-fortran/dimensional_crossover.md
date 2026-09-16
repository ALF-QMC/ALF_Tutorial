
# Exercise: Dimensional crossover

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

```{figure} ../Figures/ladder.pdf 
---
width: 60%
name: fig-ladder
---
Spin correlation functions along one leg for the Hubbard ladder. As $t_y$ grows the spin gap becomes large enough so as to detect the exponential decay of the spin correlation function on this small lattice size. The underlying physics of odd-even ladder systems is introduced in the article: Elbio Dagotto and T. M. Rice, "Surprises on the way from one- to two-dimensional quantum magnets: The ladder materials", Science **271**, 5249 (1996) pp. 618–623.
```

