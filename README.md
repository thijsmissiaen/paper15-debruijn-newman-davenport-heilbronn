# Asymptotic Shift Thresholds for Pólya–Schur Multipliers and a Certified Lower Bound on the de Bruijn–Newman Constant of the Davenport–Heilbronn Function

**DOI (Zenodo / CERN):** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24141103.svg)](https://doi.org/10.5281/zenodo.24141103) [`10.5281/zenodo.24141103`](https://doi.org/10.5281/zenodo.24141103)  
**Author:** Thijs Missiaen (Independent Researcher, Belgium · ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810))  
**Series:** Transverse Wronskian & Phase-Space Analysis of Riemann's $\Xi$-Function (Paper XV, v15.2.0)  
**MSC 2020:** Primary 11M41, 11M26; Secondary 30D15, 26C10, 65G20

---

## Overview

For the self-dual Dirichlet combination introduced by Davenport and Heilbronn (1936),

$$
\mathcal{D}_{\mathrm{DH}}(s) := \frac{1 - \mathrm{i}\kappa}{2}\thinspace L(s,\chi_5) + \frac{1 + \mathrm{i}\kappa}{2}\thinspace L(s,\overline{\chi_5}), \qquad \kappa := \frac{\sqrt{10 - 2\sqrt{5}} - 2}{\sqrt{5} - 1},
$$

the completed entire function $\Xi _ {\mathrm{DH}}(z) := \xi _ {\mathrm{DH}}(\frac{1}{2} + \frac{\mathrm{i}z}{2}) = 2\int _ 0^\infty \Phi _ {\mathrm{DH}}(u)\cos(zu)\thinspace\mathrm{d}u$ satisfies $\Xi _ {\mathrm{DH}}(z) = \Xi _ {\mathrm{DH}}(-z) = \overline{\Xi _ {\mathrm{DH}}(\bar{z})}$ yet lacks an Euler product and possesses infinitely many zeros both in $\mathrm{Re}(s) > 1$ and off the critical line inside $0 < \mathrm{Re}(s) < 1$, with off-line supremum $\sigma^{\ast} := \sup\lbrace\mathrm{Re}(\rho) : \mathcal{D} _ {\mathrm{DH}}(\rho) = 0\rbrace = 1.120362498183\dots$ (Bombieri–Ghosh 2011).

This paper establishes sufficient asymptotic shift conditions for real entire Pólya–Schur differential multipliers $G(-\partial _ z^2)$ acting on $\Xi _ {\mathrm{DH}}(z)$ and certifies explicit unconditional bounds on its de Bruijn–Newman constant $\Lambda _ {\mathrm{DH}}$:

### 1. Sufficient Stieltjes Shift Conditions and Lacunary Root-Gap Obstruction (Section 2)

Writing $G(w) = \mathrm{e}^{\lambda w}\prod _ {j=1}^\infty (1 + w/w _ j)$ ($\lambda \ge 0$, $\sum w _ j^{-1} < \infty$) and $S _ G(u) := \frac{\mathrm{d}}{\mathrm{d}u}\log G(u^2) = 4\lambda u + \sum _ {j=1}^\infty \frac{2u}{w _ j + u^2}$, we generalize Dobner's (2021) modified Dirichlet-series approximation via Pólya's conjugate-diagram Borel–Laplace representation and single-zero Diophantine translation (avoiding Righetti's 2017 non-density phenomenon):

* *Zero-shift persistence* ($\lambda = 0$, $\sum w _ j^{-1/2} < \infty$, including every polynomial and every multiplier of upper order $\rho(G) < 1/2$): $\lim _ {u\to\infty} S _ G(u) = 0$ and $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ inside $\lvert\mathrm{Im}(z)\rvert > 1$; in particular, $(1 - \beta\partial _ z^2)\Xi _ {\mathrm{DH}}(z)$ has infinitely many non-real zeros for every $\beta > 0$.
* *Finite Stieltjes shift* ($\lambda = 0$, $\lim _ {u\to\infty} S _ G(u) = \tau \in (0,\infty){}$): $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ for every $0 < \tau < 2\sigma^{\ast} - 1$, whereas $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ for every $\tau > 2\sigma^{\ast} - 1$.
* *Divergent Stieltjes shift and root-gap regularity:* Every $G$ with $\lambda > 0$ or lower order $\underline{\rho}(G) > 1/2$ satisfies $\lim _ {u\to\infty} S _ G(u) = +\infty$, and $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ holds whenever $\lambda > 0$, or whenever $\lambda = 0$ and $w _ {j+1} - w _ j = o(\sqrt{w _ j})$. Conversely, for every $\rho \in (1/2, 1)$, setting $R _ 1 := 16$ and $R _ k := 2^{2^{R _ {k-1}}}$ ($k \ge 2$), the tower-lacunary canonical product $G _ {\mathrm{lac}}(w) = \prod _ {k=1}^\infty (1 + w/R _ k^2)^{\lfloor R _ k^{2\rho}\rfloor}$ has upper order $\rho$ and lower order $0$, satisfies $\liminf _ {u\to\infty} S _ {G _ {\mathrm{lac}}}(u) = 0$, and leaves infinitely many non-real zeros in $\lvert\mathrm{Im}(z)\rvert > 1$.

### 2. Global Doubling Bound, Threshold Separation, and Conditional Critical Boundary for Pólya's Sinc Multiplier (Section 3)

Comparing $G _ {\mathrm{cosh},h}(w) = \cosh(2h\sqrt{w})$ and $G _ {\mathrm{sinc},h}(w) = \frac{\sinh(2h\sqrt{w})}{2h\sqrt{w}}$ of identical asymptotic shift $\lim _ {u\to\infty} S _ G(u) = 2h$, a singularity-free 384-bit `FLINT/Arb` 360-ball Interval–Newton continuation across the entire compact interval $h \in [0, \frac{5}{8}]$ from $\rho _ {44838}$ (plus radius $10^{-30}$ point enclosures at $h = \frac{5}{8} = 0.625$ and $h = \frac{625121}{10^6} = 0.625121$) together with a computer-free interval-splitting doubling proof ($\int _ {-h}^h = \int _ {-h}^0 + \int _ 0^h$, noting that $\lbrace G _ {\mathrm{sinc},h}\rbrace _ {h>0}$ is not a semigroup in $h$) establish independently of Section 2:

$$
h_{\mathrm{global}}^{\ast}\bigl(G_{\mathrm{cosh}}; \Xi_{\mathrm{DH}}\bigr) = \sigma^{\ast} - \frac{1}{2} < \frac{5}{8} < h_{\mathrm{global}}^{\ast}\bigl(G_{\mathrm{sinc}}; \Xi_{\mathrm{DH}}\bigr) \le 2\sigma^{\ast} - 1 = 1.240724996\dots
$$

At the critical asymptotic boundary $h _ {\ast} = \sigma^{\ast} - \frac{1}{2}$ ($\tau = 2\sigma^{\ast} - 1$), Proposition 3.3 proves conditionally under RH for $\zeta(s)$ and GRH for $L(s, \chi _ {5,\mathrm{quad}})$ that $8 \le N _ {\mathrm{off}}(G _ {\mathrm{sinc},h _ {\ast}}(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ (with $N _ {\mathrm{off}} \ge 8$ certified unconditionally via $\rho _ {2442}$ and $\rho _ {44838}$).

### 3. Certified 12-Digit Rational Lower Bound and Ki–Kim–Lee Upper Bound on $\Lambda _ {\mathrm{DH}}$ (Sections 4–5)

Adapting the large-height asymptotic simplicity and finite-box hyperbola method of Ki, Kim, and Lee (2009) to $H _ t^{\mathrm{DH}}(z)$ establishes the strict upper bound $\Lambda _ {\mathrm{DH}} < 2(\sigma^{\ast} - \frac{1}{2})^2$ (contrasting Dobner's Riemann–Siegel-truncated modified Dirichlet polynomial $\ell^1$ barrier in Remark 4.5), while a 512-bit `FLINT/Arb` Interval–Newton enclosure with third-order derivative Taylor bound at exact rational $t = \frac{28825920237}{125000000000} = 0.230607361896$ near $\rho _ {240}$ ($t _ {\mathrm{col}}(\rho _ {240}) \approx 0.2306073618964546$) proves unconditionally:

$$
\frac{28825920237}{125000000000} < \Lambda_{\mathrm{DH}} < 2\left(\sigma^{\ast} - \frac{1}{2}\right)^2 = 0.769699258304\dots
$$

(For context, Platt and Trudgian proved $\Lambda _ \zeta \le \frac{1}{5} = 0.2 < \frac{28825920237}{125000000000}$ under the same heat-equation normalization.)

---

## Repository Structure & Reproducibility

* [`paper/paper15_debruijn_newman_dh.tex`](paper/paper15_debruijn_newman_dh.tex) – AMS-LaTeX manuscript.
* [`paper/paper15_debruijn_newman_dh.pdf`](paper/paper15_debruijn_newman_dh.pdf) – Compiled PDF (20 pages, AMS-LaTeX / pdfTeX).
* [`scripts/verify_paper15_debruijn_newman_dh.py`](scripts/verify_paper15_debruijn_newman_dh.py) – Rigorous 320-bit, 384-bit, and 512-bit `python-flint` (`FLINT/Arb`) ball-arithmetic verification suite (pass `--regenerate-figures` to rebuild `figures/`).
* [`figures/`](figures/) – Vector PDF and 300-DPI PNG figures (`fig1_heat_flow_separation_witnesses`, `fig2_polya_schur_sinc_vs_cosh`).
* [`.agents/skills/debruijn-newman-davenport-heilbronn/SKILL.md`](.agents/skills/debruijn-newman-davenport-heilbronn/SKILL.md) – Self-contained mathematical knowledge skill.

```bash
python scripts/verify_paper15_debruijn_newman_dh.py
```
