# Monte-Carlo pitfalls and comparison to exact results

This chapter is devoted to a discussion of Monte Carlo pitfalls and how to trust your error-bars.  For each new implementation of a model Hamiltonian it is imperative to verify  the result. This can only be done by comparing to exact results on small lattices.  This allows testing for all possible sources systematic errors, including  bad random number generators.
Here, we will first explore the effects of the systematic Trotter error using the Hubbard model on a 2D square lattice. Then, we will compare results from ALF to exact diagonalisation via the four-site Hubbard chain.
Eventually, we consider the warm-up time and autocorrelation and its effects on results.

 
