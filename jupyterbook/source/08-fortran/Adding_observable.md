
# Exercise: Adding a new observable

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

In [solution](Solutions/Exercise_3/Hamiltonian_Hubbard_Plain_Vanilla_smod.F90-Exercise_3) we have the modified submodule `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90`, as well as the [parameter](Solutions/Exercise_3/Start/parameters) file. The following are the essential steps to be carried out:

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

