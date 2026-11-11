# Math notes

## Modified plasma dispersion function

For a kappa distribution along ``B_0`` with index ``\kappa_\|``, drift ``u_0`` and width ``\theta_\|``, the Landau integral becomes the modified plasma dispersion function of Summers and Thorne,

```math
Z_{\kappa_\|}(\xi_n)=\frac{\Gamma(\kappa_\|+1)}{\sqrt{\pi \kappa_\|^3}\,\Gamma(\kappa_\|-1/2)} \int_{-\infty}^{\infty} \frac{dx}{(x-\xi_n)\,(1+x^2/\kappa_\|)^{\kappa_\|+1}}, \qquad \Im\xi_n>0,
```

with ``x=(v_\|-u_0)/\theta_\|`` and ``\xi_n=(\omega-n\Omega-k_\| u_0)/(k_\|\theta_\|)``, where ``\Omega`` is the signed gyrofrequency and ``n`` the cyclotron harmonic. For integer ``\kappa_\|`` the integrand is rational, and the integral closes to

```math
Z_{\kappa_\|}(\xi_n)=\frac{i(\kappa_\|-1/2)}{2\kappa_\|^{3/2}}\,\frac{\kappa_\|!}{(2\kappa_\|)!} \sum_{l=1}^{\kappa_\|+1} \frac{(2\kappa_\|-l+1)!}{(\kappa_\|-l+1)!}\left(\frac{2i\sqrt{\kappa_\|}}{\xi_n+i\sqrt{\kappa_\|}}\right)^{l},
```

valid for all ``\xi_n \neq -i\sqrt{\kappa_\|}``, which also continues it analytically to ``\Im\xi_n \le 0``. In terms of ``\omega``, ``\xi_n + i\sqrt{\kappa_\|} = (\omega - c_n)/(k_\|\theta_\|)`` with

```math
c_n = n\Omega + k_\|\left(u_0 - i\sqrt{\kappa_\|}\,\theta_\|\right),
```

so the ``l``-th term is ``\propto (2i\sqrt{\kappa_\|}\,\theta_\| k_\|)^l/(\omega-c_n)^l``. ``Z_{\kappa_\|}`` is therefore exactly a single pole of order ``\kappa_\|+1`` at ``\omega = c_n``, lying ``\sqrt{\kappa_\|}\,\theta_\| k_\|`` below the real axis. This holds for ``k_\| > 0``; `BOPBK` maps ``k_\| < 0`` onto it by mirroring ``v_\| \to -v_\|`` and reversing the drifts.

## Linear eigenvalue problem (`BOPBK`)

`BOPBK` turns the dispersion relation into an eigenvalue problem ``\omega X = M X`` whose eigenvalues are the roots ``\omega``. Units have ``c = \varepsilon_0 = 1``, and quantities are per species unless summed.

### The current as a sum of poles

Substituting the closed form into the conductivity, harmonic ``n`` contributes to the current only through powers of ``(\omega - c_n)^{-1}``. Each term ``l`` of ``Z_{\kappa_\|}`` enters twice: as ``(\omega-c_n)^{-l}`` with a 3×3 coefficient ``H_{n,l}`` (`_pbk_head`), and one order higher, as ``(\omega-c_n)^{-(l+1)}``, with ``T_{n,l}`` (`_pbk_tail`). Collecting equal powers, with ``m = \kappa_\|+2`` the highest order reached,

```math
\omega J = i\omega_p^2 E + \sum_{n=-N}^{N} \sum_{j=1}^{m} \frac{P_{n,j}\,E}{(\omega - c_n)^j},
\qquad P_{n,j} = H_{n,j} + T_{n,j-1},\quad H_{n,m} = T_{n,0} = 0 .
```

Both coefficients of term ``l`` carry ``(2i\sqrt{\kappa_\|}\,\theta_\| k_\|)^{l-1}`` (`_csl`), so ``\lVert P_{n,j}\rVert \propto s^{j-1}`` with the pole width

```math
s = 2\sqrt{\kappa_\|}\,\theta_\|\,\lvert k_\|\rvert = 2\lvert\Im c_n\rvert .
```

### Unfolding each pole into a chain

A pole of order ``m`` is not linear in ``\omega``, but it can be peeled one order at a time. For each harmonic, species and field component, introduce auxiliary states

```math
y_j = \sum_{i=j}^{m} \frac{P_{n,i}\,E}{(\omega - c_n)^{i-j+1}}, \qquad j = 1,\dots,m :
```

``y_j`` holds orders ``j`` through ``m`` of the expansion, each lowered by ``j-1``. So ``y_1`` is exactly the harmonic's term in the current, and successive states differ by one factor of the pole, ``(\omega - c_n)\,y_j = y_{j+1} + P_{n,j}E``. Written out, this is linear in ``\omega``:

```math
\omega\, y_j = c_n y_j + y_{j+1} + P_{n,j} E, \qquad y_{m+1} = 0 .
```

Each chain is thus a Jordan block at ``c_n``, driven by ``E`` and read out through ``y_1``. One chain holds every term of the expansion because all terms share the pole ``c_n``.

### The full system

The chains, the current, and Maxwell's curl equations close the system:

```math
\begin{aligned}
\omega\, y_j &= c_n y_j + y_{j+1} + P_{n,j} E, \\
\omega J &= i\omega_p^2 E + \textstyle\sum y_1, \\
\omega E &= -\mathbf{k}\times B - iJ, \\
\omega B &= \mathbf{k}\times E,
\end{aligned}
```

where ``\sum y_1`` runs over species and harmonics, and ``\omega_p^2`` is summed over species. Eliminating ``y``, ``J`` and ``B`` recovers the dispersion relation ``\det D(\omega,\mathbf{k}) = 0``, so the eigenvalues of ``M`` give all its roots at once.

### Scaling the chain to the pole width

As written, the chain is badly conditioned. `BOPBK` stores ``x_j = s^{-(j-1)} y_j`` instead:

```math
\omega\, x_j = c_n x_j + s\, x_{j+1} + s^{-(j-1)} P_{n,j} E .
```

This is the diagonal similarity ``D^{-1}MD`` with ``D = \operatorname{diag}(s^{j-1})``, and ``x_1 = y_1``, so the eigenvalues and the current are unchanged. What changes is their sensitivity to roundoff.

A Jordan block of size ``m`` with eigenvalue ``c_n`` and superdiagonal ``a``, perturbed by ``\delta`` in its lower-left corner, has ``(\mu - c_n)^m = \delta\, a^{m-1}``. Its eigenvalue spreads to radius

```math
\lvert \mu - c_n \rvert = a\,(\delta/a)^{1/m},
```

and eigenvalues of the full matrix near ``c_n`` inherit this sensitivity. The scale that matters is ``s``: roots near the pole lie at distance ``\sim s`` from it, and the pole lies ``s/2`` below the real axis. Relative to ``s``, the spread is

- unscaled, ``a = 1``: ``\delta^{1/m}/s``, unbounded as ``s \to 0``. The couplings ``P_{n,j} \propto s^{j-1}`` also span many orders of magnitude, which LAPACK's power-of-2 balancing does not detect. Once the spread exceeds ``s/2``, roots cross the real axis and show up as spurious growth, ``\Im\omega > 0``.
- scaled, ``a = s``: ``(\delta/s)^{1/m}``, with every coupling ``s^{-(j-1)}P_{n,j} = O(1)``. For roundoff ``\delta \sim \varepsilon s``, with ``\varepsilon`` the machine epsilon, the spread is ``\varepsilon^{1/m} s``, independent of ``k_\|``.

The remaining factor ``\varepsilon^{1/m}`` still grows with ``\kappa_\|``: at ``m = 26`` it is ``\approx 0.25``, so the spread is already half the pole depth ``s/2``. `BOPBK` therefore rejects ``\kappa_\| > 24``.

At ``k_\| = 0`` the pole width vanishes and `BOPBK` keeps ``a = 1``. The matrix jumps there, but its eigenvalues stay continuous.
