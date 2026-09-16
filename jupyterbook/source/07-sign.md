# My very first sign problem

In this chapter, we will use a notebook to compute the average sign of the
square-lattice Hubbard model at $U/t=4$ as a function of the chemical potential
and lattice size. We will see that simulations can remain feasible in the
presence of a sign problem, provided that the average sign is not too small.
We will also explore how strongly the sign problem depends on the choice of
Hubbard--Stratonovich transformation.

## Complex configuration weights and reweighting

After introducing the auxiliary fields and integrating out the fermions, the
partition function can be written as

$$
Z
=
\sum_C W_{\mathrm{B}}(C)W_{\mathrm{F}}(C),
$$

where $C$ denotes an auxiliary-field configuration. The factor
$W_{\mathrm{B}}(C)$ contains the purely bosonic contribution, while
$W_{\mathrm{F}}(C)$ results from tracing out the fermionic degrees of freedom
and is generally expressed in terms of fermionic determinants.

For the following discussion, we assume that the bosonic weight is
non-negative,

$$
W_{\mathrm{B}}(C)\geq0.
$$

In general, however, the fermionic weight can be complex. The complete
configuration weight can then no longer be interpreted as a probability
distribution. We can nevertheless formulate a reweighting scheme by writing
the fermionic weight in terms of its absolute value and phase,

$$
W_{\mathrm{F}}(C)
=
\left|W_{\mathrm{F}}(C)\right|e^{i\theta(C)}.
$$

The partition function consequently takes the form

$$
Z
=
\sum_C
W_{\mathrm{B}}(C)
\left|W_{\mathrm{F}}(C)\right|
e^{i\theta(C)}.
$$

To construct a valid sampling probability, we define the absolute-weight
partition function

$$
Z_{\mathrm{abs}}
\equiv
\sum_C
W_{\mathrm{B}}(C)
\left|W_{\mathrm{F}}(C)\right|.
$$

In contrast to the physical partition function, every contribution to
$Z_{\mathrm{abs}}$ is real and non-negative. The corresponding normalized
probability distribution is

$$
P_{\mathrm{abs}}(C)
\equiv
\frac{
W_{\mathrm{B}}(C)\left|W_{\mathrm{F}}(C)\right|
}{
Z_{\mathrm{abs}}
}.
$$

Using this probability distribution, the physical partition function can be
rewritten as

$$
\begin{aligned}
Z
&=
\sum_C
W_{\mathrm{B}}(C)
\left|W_{\mathrm{F}}(C)\right|
e^{i\theta(C)}
\\
&=
Z_{\mathrm{abs}}
\sum_C
P_{\mathrm{abs}}(C)e^{i\theta(C)}
\\
&=
Z_{\mathrm{abs}}
\left\langle e^{i\theta(C)}\right\rangle_{\mathrm{abs}}.
\end{aligned}
$$

Here,

$$
\left\langle e^{i\theta(C)}\right\rangle_{\mathrm{abs}}
\equiv
\sum_C
P_{\mathrm{abs}}(C)e^{i\theta(C)}
$$

is the average phase.

The thermal expectation value of an observable $\hat{O}$ is

$$
\langle\hat{O}\rangle
=
\frac{
\displaystyle
\sum_C
W_{\mathrm{B}}(C)W_{\mathrm{F}}(C)
\langle\!\langle\hat{O}\rangle\!\rangle_C
}{
\displaystyle
\sum_C
W_{\mathrm{B}}(C)W_{\mathrm{F}}(C)
},
$$

where $\langle\!\langle\hat{O}\rangle\!\rangle_C$ denotes the value of the
observable for a fixed auxiliary-field configuration. Using
$W_{\mathrm{F}}(C)=|W_{\mathrm{F}}(C)|e^{i\theta(C)}$, the numerator becomes

$$
\begin{aligned}
&\sum_C
W_{\mathrm{B}}(C)W_{\mathrm{F}}(C)
\langle\!\langle\hat{O}\rangle\!\rangle_C
\\
&\qquad =
Z_{\mathrm{abs}}
\sum_C
P_{\mathrm{abs}}(C)e^{i\theta(C)}
\langle\!\langle\hat{O}\rangle\!\rangle_C
\\
&\qquad =
Z_{\mathrm{abs}}
\left\langle
e^{i\theta(C)}
\langle\!\langle\hat{O}\rangle\!\rangle_C
\right\rangle_{\mathrm{abs}}.
\end{aligned}
$$

Together with

$$
Z
=
Z_{\mathrm{abs}}
\left\langle e^{i\theta(C)}\right\rangle_{\mathrm{abs}},
$$

this yields the reweighting formula

$$
\langle\hat{O}\rangle
=
\frac{
\left\langle
e^{i\theta(C)}
\langle\!\langle\hat{O}\rangle\!\rangle_C
\right\rangle_{\mathrm{abs}}
}{
\left\langle e^{i\theta(C)}\right\rangle_{\mathrm{abs}}
}.
$$

The equations above describe the general phase problem. In the simulations
considered below, the fermionic weight is real. Its phase is therefore
restricted to

$$
\begin{cases}
\theta(C)\equiv0\pmod{2\pi},
& W_{\mathrm{F}}(C)>0,\\
\theta(C)\equiv\pi\pmod{2\pi},
& W_{\mathrm{F}}(C)<0.
\end{cases}
$$

We can consequently define the sign

$$
s(C)
\equiv
e^{i\theta(C)}
=
\operatorname{sgn}\!\left[W_{\mathrm{F}}(C)\right]
\in\{-1,+1\}.
$$

The average phase then reduces to the average sign,

$$
\left\langle e^{i\theta(C)}\right\rangle_{\mathrm{abs}}
=
\left\langle s(C)\right\rangle_{\mathrm{abs}},
$$

and the reweighting formula becomes

$$
\langle\hat{O}\rangle
=
\frac{
\left\langle
s(C)\langle\!\langle\hat{O}\rangle\!\rangle_C
\right\rangle_{\mathrm{abs}}
}{
\left\langle s(C)\right\rangle_{\mathrm{abs}}
}.
$$

## Why is the sign problem difficult?

At first sight, it seems that we have solved the problem by sampling with the
positive distribution $P_{\mathrm{abs}}(C)$. However, the sign has not
disappeared. It enters both the numerator and the denominator of the
reweighted expectation value. If the average sign becomes very small, these
quantities become increasingly difficult to estimate accurately. Reweighting
therefore allows us to perform the sampling, but it does not remove the sign
problem.

To see this more explicitly, consider the relative statistical uncertainty of
the sampled average sign. Since $s(C)^2=1$, it is given by

$$
\frac{
\Delta\left\langle s\right\rangle_{\mathrm{abs}}
}{
\left|\left\langle s\right\rangle_{\mathrm{abs}}\right|
}
=
\frac{
\sqrt{1-\left\langle s\right\rangle_{\mathrm{abs}}^2}
}{
\sqrt{N}
\left|\left\langle s\right\rangle_{\mathrm{abs}}\right|
},
$$

where $N$ denotes the number of statistically independent configurations. The
average sign generally decreases exponentially with inverse temperature and
system size,

$$
\left|\left\langle s\right\rangle_{\mathrm{abs}}\right|
\sim
e^{-\beta N_s\Delta f},
$$

where $N_s$ is the number of lattice sites and $\Delta f\geq0$ is the
difference between the free-energy densities of the physical and
absolute-weight ensembles. For a small average sign, the relative uncertainty
therefore scales as

$$
\frac{
\Delta\left\langle s\right\rangle_{\mathrm{abs}}
}{
\left|\left\langle s\right\rangle_{\mathrm{abs}}\right|
}
\sim
\frac{e^{\beta N_s\Delta f}}{\sqrt{N}}.
$$

Maintaining a fixed relative uncertainty consequently requires

$$
N
\sim
e^{2\beta N_s\Delta f}.
$$

The required computational effort therefore grows exponentially with inverse
temperature and system size. Reweighting is exact, but the cancellations
between positive and negative configurations make the average sign
exponentially difficult to resolve. The same cancellation mechanism applies
to complex weights, with the average sign replaced by the average phase.



**TODO and Planned numerical part**

- pyALF: Hubbard model at $U/t=4$ for different values of $\beta$, lattice
  sizes and HS transformations.
- Plot of the average sign for the students to reproduce.
- Project idea: overdoped Hubbard model in metallic regime with a manageable sign problem;
  examine the equal-time Green function near the Fermi surface and the momentum
  occupation.


:::{figure} Figures/Average_sign_example_plot.jpeg
:width: 600px
:align: center
:alt: Average sign of the Hubbard model

Average sign as a function of the chemical potential for different system sizes.
:::