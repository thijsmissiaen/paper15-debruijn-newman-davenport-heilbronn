# Stieltjes Regimes of Pólya–Schur Multipliers and Strict de Bruijn–Newman Bounds for the Davenport–Heilbronn Function

**Author:** Thijs Missiaen (Independent Researcher, Belgium · ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810))  
**Series:** Transverse Wronskian & Phase-Space Analysis of Riemann's $\Xi$-Function (Paper XV)  
**MSC 2020:** Primary 11M41, 11M26; Secondary 30D15, 26C10, 65G20

---

## Overview

For the self-dual Dirichlet combination introduced by Davenport and Heilbronn (1936),

$$\mathcal{D}_{\mathrm{DH}}(s) := \frac{1 - \mathrm{i}\kappa}{2}\,L(s,\chi_5) + \frac{1 + \mathrm{i}\kappa}{2}\,L(s,\overline{\chi_5}), \qquad \kappa := \frac{\sqrt{10 - 2\sqrt{5}} - 2}{\sqrt{5} - 1},$$

the completed entire function $\Xi _ {\mathrm{DH}}(z) := \xi _ {\mathrm{DH}}(\frac{1}{2} + \frac{\mathrm{i}z}{2}) = 2\int _ 0^\infty \Phi _ {\mathrm{DH}}(u)\cos(zu)\,\mathrm{d}u$ satisfies the real-even functional equation $\Xi _ {\mathrm{DH}}(z) = \Xi _ {\mathrm{DH}}(-z) = \overline{\Xi _ {\mathrm{DH}}(\bar{z})}$ yet lacks an Euler product and possesses infinitely many zeros both in $\mathrm{Re}(s) > 1$ and off the critical line inside $0 < \mathrm{Re}(s) < 1$, with Bohr–Kronecker supremum $\sigma^* := \sup\{\mathrm{Re}\,\rho : \mathcal{D} _ {\mathrm{DH}}(\rho) = 0\} = 1.120362498183\dots$ (Bombieri–Ghosh 2011).

This paper classifies the action of real entire Pólya–Schur differential multipliers $G(-\partial _ z^2)$ on $\Xi _ {\mathrm{DH}}(z)$ and establishes strict unconditional two-sided bounds on its de Bruijn–Newman constant $\Lambda _ {\mathrm{DH}}$:

1. **Three-Regime Stieltjes Classification and Root-Gap Obstruction (Section 2):** Writing $G(w) = \mathrm{e}^{\lambda w}\prod _ {j=1}^\infty (1 + w/w _ j)$ ($\lambda \ge 0$, $\sum w _ j^{-1} < \infty$) and $S _ G(u) := \frac{\mathrm{d}}{\mathrm{d}u}\log G(u^2) = 4\lambda u + \sum _ {j=1}^\infty \frac{2u}{w _ j + u^2}$, we prove via Pólya's conjugate-diagram Borel–Laplace representation and Rolle–Hurwitz subproduct monotonicity:
   * *Regime I* ($\lambda = 0$, $\sum w _ j^{-1/2} < \infty$, including every polynomial and every multiplier of upper order $\rho(G) < 1/2$): $\lim _ {u\to\infty} S _ G(u) = 0$ and $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ inside $\lvert\mathrm{Im}\,z\rvert > 1$; in particular, $(1 - \beta\partial _ z^2)\Xi _ {\mathrm{DH}}(z)$ has infinitely many non-real zeros for every $\beta > 0$.
   * *Regime II* ($\lambda = 0$, $\lim _ {u\to\infty} S _ G(u) = \tau \in (0,\infty){}$): $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ for every $0 < \tau < 2\sigma^* - 1$, whereas $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ for every $\tau > 2\sigma^* - 1$.
   * *Regime III:* Every $G$ with $\lambda > 0$ or lower order $\underline{\rho}(G) > 1/2$ satisfies $\lim _ {u\to\infty} S _ G(u) = +\infty$, and $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ holds whenever $\lambda > 0$, or whenever $\lambda = 0$ and $w _ {j+1} - w _ j = o(\sqrt{w _ j})$. Conversely, for every $\rho \in (1/2, 1)$, setting $R _ 1 := 16$ and $R _ k := 2^{2^{R _ {k-1}}}$ ($k \ge 2$), the tower-lacunary canonical product $G _ {\mathrm{lac}}(w) = \prod _ {k=1}^\infty (1 + w/R _ k^2)^{\lfloor R _ k^{2\rho}\rfloor}$ has upper order $\rho$ and lower order $0$, satisfies $\liminf _ {u\to\infty} S _ {G _ {\mathrm{lac}}}(u) = 0$, and leaves infinitely many non-real zeros in $\lvert\mathrm{Im}\,z\rvert > 1$.
2. **Global versus Asymptotic Threshold Separation for Pólya's Sinc Multiplier (Section 3):** Comparing $G _ {\mathrm{cosh},h}(w) = \cosh(2h\sqrt{w})$ and $G _ {\mathrm{sinc},h}(w) = \frac{\sinh(2h\sqrt{w})}{2h\sqrt{w}}$ of identical asymptotic shift $\lim _ {u\to\infty} S _ G(u) = 2h$, a 384-bit `FLINT/Arb` Taylor-jet + 512-arc Cauchy enclosure at $h = \frac{311}{500}$ together with a 64-ball parametric Interval–Newton cover across $[\frac{310}{500}, \frac{311}{500}] \supset [\sigma^* - \frac{1}{2}, \frac{311}{500}]$ near $\rho _ {2442}$ proves unconditionally without any semigroup assumption:

$$h_{\mathrm{global}}^*\bigl(G_{\mathrm{cosh}}; \Xi_{\mathrm{DH}}\bigr) = \sigma^* - \frac{1}{2} < \frac{311}{500} < h_{\mathrm{global}}^*\bigl(G_{\mathrm{sinc}}; \Xi_{\mathrm{DH}}\bigr).$$

3. **Shifted Gauss–Weierstrass Asymptotics and Strict Two-Sided Bounds on $\Lambda _ {\mathrm{DH}}$ (Sections 4–5):** Folding $H _ t^{\mathrm{DH}}(2\gamma - 2\mathrm{i}\delta) = A _ t(\gamma,\delta) + \overline{A _ t(\gamma,-\delta)}$ on $[0,\infty)$, shifting $v \mapsto u + \mathrm{i}\frac{\pi t}{8}$ to cancel Stirling's linear phase identically, and applying de Bruijn's equilateral-hyperbola contraction on $\lvert\mathrm{Re}\,z\rvert \le 2\gamma _ 0(t _ 1)$ proves $\Lambda _ {\mathrm{DH}} < 2(\sigma^* - \frac{1}{2})^2$. Combining 640-bit, 704-bit, and 1350-bit `FLINT/Arb` Interval–Newton off-real zero certificates at $t = \frac{1}{5}$ ($\rho _ {240}$ and $\rho _ {531}$) and at $t = \frac{1153}{5000} = 0.2306$ ($\rho _ {240}$) with Platt and Trudgian's bound $\Lambda _ \zeta \le \frac{1}{5}$ establishes unconditionally:

$$\Lambda_\zeta \le \frac{1}{5} < \frac{1153}{5000} < \Lambda_{\mathrm{DH}} < 2\left(\sigma^* - \frac{1}{2}\right)^2 = 0.769699258304\dots$$

---

## Repository Structure & Reproducibility

* [`paper/paper15_debruijn_newman_dh.tex`](paper/paper15_debruijn_newman_dh.tex) – AMS-LaTeX manuscript.
* [`paper/paper15_debruijn_newman_dh.pdf`](paper/paper15_debruijn_newman_dh.pdf) – Compiled PDF (14 pages, AMS-LaTeX / pdfTeX).
* [`scripts/verify_paper15_debruijn_newman_dh.py`](scripts/verify_paper15_debruijn_newman_dh.py) – Rigorous 320-bit to 1350-bit `python-flint` (`FLINT/Arb`) ball-arithmetic verification and vector figure generation suite.
* [`figures/`](figures/) – Vector PDF and 300-DPI PNG figures (`fig1_heat_flow_separation_witnesses`, `fig2_polya_schur_sinc_vs_cosh`).
* [`.agents/skills/debruijn-newman-davenport-heilbronn/SKILL.md`](.agents/skills/debruijn-newman-davenport-heilbronn/SKILL.md) – Self-contained mathematical knowledge skill.

```bash
python scripts/verify_paper15_debruijn_newman_dh.py
```
