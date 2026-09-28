# What is the sign problem?

## Configuration weights and reweighting

After introducing the auxiliary fields and integrating out the fermions, the partition function can be written as

$$
Z
=
\sum_C e^{-S(C)}
\equiv
\sum_C \mathcal{W}(C),
$$

where $C$ denotes an auxiliary-field configuration, $S(C)$ is the corresponding effective action, and

$$
\mathcal{W}(C) \equiv e^{-S(C)}
$$

is its configuration weight. In general, the effective action and hence the configuration weight can be complex. The weights can therefore no longer be interpreted as probabilities, since a probability distribution must be real and non-negative. Consequently, they cannot be used directly for importance sampling.

To construct a suitable sampling probability, we follow the reweighting scheme
used in the ALF documentation. Since the partition function is real, we can
write

$$
Z
=
\operatorname{Re}Z
=
\sum_C \operatorname{Re}\!\left[e^{-S(C)}\right].
$$

Although $\operatorname{Re}[e^{-S(C)}]$ is real, it is not necessarily
positive. We therefore define the sign of a configuration as

$$
\operatorname{sgn}(C)
\equiv
\frac{
\operatorname{Re}\!\left[e^{-S(C)}\right]
}{
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
}
\in \{-1,+1\}.
$$

A non-negative probability distribution can now be constructed from the
absolute value of the real part of the configuration weight,

$$
\bar{P}(C)
\equiv
\frac{
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
}{
\displaystyle
\sum_C
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
}.
$$

For any configuration-dependent quantity $X(C)$, we denote averages with
respect to this probability distribution by

$$
\left\langle X\right\rangle_{\bar{P}}
\equiv
\sum_C \bar{P}(C)X(C).
$$

The average sign is therefore given by

$$
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
=
\frac{
\displaystyle
\sum_C
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
\operatorname{sgn}(C)
}{
\displaystyle
\sum_C
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
}.
$$

Let $\langle\!\langle\hat{O}\rangle\!\rangle_C$ denote the estimator of an
observable $\hat{O}$ for a fixed auxiliary-field configuration. Its thermal
expectation value is

$$
\left\langle\hat{O}\right\rangle
=
\frac{
\displaystyle
\sum_C
e^{-S(C)}
\langle\!\langle\hat{O}\rangle\!\rangle_C
}{
\displaystyle
\sum_C e^{-S(C)}
}.
$$

Using

$$
\operatorname{Re}\!\left[e^{-S(C)}\right]
=
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
\operatorname{sgn}(C),
$$

we can rewrite this expression as

$$
\begin{aligned}
\left\langle\hat{O}\right\rangle
&=
\frac{
\displaystyle
\sum_C
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
\operatorname{sgn}(C)
\frac{e^{-S(C)}}{\operatorname{Re}[e^{-S(C)}]}
\langle\!\langle\hat{O}\rangle\!\rangle_C
}{
\displaystyle
\sum_C
\left|\operatorname{Re}\!\left[e^{-S(C)}\right]\right|
\operatorname{sgn}(C)
}
\\[1ex]
&=
\frac{
\displaystyle
\left\langle
\operatorname{sgn}(C)
\frac{e^{-S(C)}}{\operatorname{Re}[e^{-S(C)}]}
\langle\!\langle\hat{O}\rangle\!\rangle_C
\right\rangle_{\bar{P}}
}{
\displaystyle
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
}.
\end{aligned}
$$

The factor

$$
\frac{e^{-S(C)}}{\operatorname{Re}[e^{-S(C)}]}
$$

ensures that the information contained in the original complex configuration
weight is retained. If $e^{-S(C)}$ is real for every configuration, this factor
equals one and the expression reduces to the usual sign-reweighting formula.

By introducing the probability distribution $\bar{P}(C)$, we can sample configurations using a real and non-negative weight. This makes Monte Carlo sampling possible, but it does not eliminate the sign
problem.

## Why is the sign problem difficult?

The sign still enters both the numerator and the denominator of the reweighted
expectation value. If the average sign becomes very small, both quantities
become increasingly difficult to estimate accurately.

To see this more explicitly, consider the relative statistical uncertainty of
the sampled average sign. Since

$$
\operatorname{sgn}(C)^2=1,
$$

it is given by

$$
\frac{
\Delta\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
}{
\left|
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
\right|
}
=
\frac{
\sqrt{
1-
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}^2
}
}{
\sqrt{N}
\left|
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
\right|
},
$$

where $N$ denotes the number of statistically independent configurations.

The normalization of the sampling distribution,

$$
Z_{\bar{P}}
\equiv
\sum_C
\left|
\operatorname{Re}\!\left[e^{-S(C)}\right]
\right|,
$$

can be interpreted as an auxiliary partition function. The average sign is
the ratio of the physical and auxiliary partition functions,

$$
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
=
\frac{Z}{Z_{\bar{P}}}.
$$

The average sign generally becomes exponentially small as the inverse
temperature and system size increase [@10.1103/PhysRevB.41.9301],

$$
\left|
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
\right|
\sim
e^{-\Delta f\beta N_s},
$$

where $N_s$ is the number of lattice sites and

$$
\Delta f \equiv f-f_{\bar{P}}\geq0
$$

is the difference between the free-energy density $f$ of the physical system
and the free-energy density $f_{\bar{P}}$ associated with the auxiliary
partition function $Z_{\bar{P}}$. Its value depends on the model and the chosen
formulation.

For a small average sign, the relative uncertainty therefore scales as

$$
\frac{
\Delta\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
}{
\left|
\left\langle\operatorname{sgn}\right\rangle_{\bar{P}}
\right|
}
\sim
\frac{e^{\Delta f \beta N_s}}{\sqrt{N}}.
$$

Maintaining a fixed relative uncertainty consequently requires

$$
N
\sim
e^{2\Delta f \beta N_s}.
$$

The required computational effort therefore grows exponentially with inverse
temperature and system size. The reweighting procedure itself introduces no
additional approximation, but cancellations between configurations with
positive and negative real weights make the average sign exponentially
difficult to resolve. More generally, finding a generic solution to the
fermion sign problem has been shown to be NP-hard
[@10.1103/PhysRevLett.94.170201].

Next, we consider the square-lattice Hubbard model and investigate how the sign
problem develops away from half filling at finite chemical potential.