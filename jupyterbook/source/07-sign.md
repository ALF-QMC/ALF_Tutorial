# My very first sign problem

In this chapter, we use a notebook to compute the average sign of the
square-lattice Hubbard model at $U/t=4$ as a function of chemical potential,
inverse temperature, and lattice size. We will see that simulations can remain feasible in the
presence of a sign problem, provided that the average sign is not too small.
We will also explore how strongly the sign problem depends on the choice of
Hubbard-Stratonovich transformation.

The first section introduces complex configuration weights, reweighting, and
the exponential scaling associated with a small average sign. The second
section uses pyALF to investigate the average sign numerically.