---
name: debruijn-newman-davenport-heilbronn
description: Canonical definitions, Stieltjes Polya-Schur multiplier regimes, FLINT/Arb interval ball certificates, and scope boundaries for Paper XV on the de Bruijn-Newman constant and Polya-Schur multipliers for the Davenport-Heilbronn function.
---

# Skill: Pólya–Schur Regimes and Strict $\Lambda _ {\mathrm{DH}}$ Bounds (Paper XV)

## 1. Canonical Normalization & Functional Equation Folding

For the Davenport–Heilbronn combination $\mathcal{D} _ {\mathrm{DH}}(s) = \frac{1 - \mathrm{i}\kappa}{2}L(s,\chi _ 5) + \frac{1 + \mathrm{i}\kappa}{2}L(s,\overline{\chi _ 5})$ with $\kappa = \frac{\sqrt{10 - 2\sqrt{5}} - 2}{\sqrt{5} - 1}$, the completed function under the canonical de Bruijn–Newman normalization is:

$$
\xi_{\mathrm{DH}}(s) := \left(\frac{5}{\pi}\right)^{(s+1)/2}\Gamma\left(\frac{s+1}{2}\right)\mathcal{D}_{\mathrm{DH}}(s), \qquad \Xi_{\mathrm{DH}}(z) := \xi_{\mathrm{DH}}\left(\frac{1}{2} + \frac{\mathrm{i}z}{2}\right) = 2\int_0^\infty \Phi_{\mathrm{DH}}(u)\cos(zu)\thinspace\mathrm{d}u.
$$

Under backward heat evolution $H _ t^{\mathrm{DH}}(z) = \mathrm{e}^{-t\partial _ z^2}\Xi _ {\mathrm{DH}}(z) = 2\int _ 0^\infty \mathrm{e}^{tu^2}\Phi _ {\mathrm{DH}}(u)\cos(zu)\thinspace\mathrm{d}u$, writing $z = 2\gamma - 2\mathrm{i}\delta$ ($\gamma,\delta \in \mathbb{R}$) and folding $\xi _ {\mathrm{DH}}(\frac{1}{2}-v+\mathrm{i}\gamma) = \overline{\xi _ {\mathrm{DH}}(\frac{1}{2}+v+\mathrm{i}\gamma)}$ yields the exact half-line Gauss–Weierstrass representation and boundary identity:

$$
H_t^{\mathrm{DH}}(2\gamma - 2\mathrm{i}\delta) = A_t(\gamma,\delta) + \overline{A_t(\gamma,-\delta)}, \qquad A_t(\gamma,\delta) := \frac{1}{\sqrt{\pi t}}\int_0^\infty \mathrm{e}^{-(v-\delta)^2/t}\thinspace\xi_{\mathrm{DH}}\left(\frac{1}{2} + v + \mathrm{i}\gamma\right)\mathrm{d}v.
$$

## 2. Certified Unconditional Results

### 2.1 Bohr–Kronecker Supremum Enclosure (Lemma 3.1)

Solving Bombieri and Ghosh's prime-angle equation via 320-bit Arb Möbius inversion and the Intermediate Value Theorem gives $\sigma^{\ast} = 1.12036249818332506180145892226464\dots$, $h _ {\infty}^{\ast} = \sigma^{\ast} - \frac{1}{2} = 0.62036249818332506180145892226464\dots$, and $2(\sigma^{\ast} - \frac{1}{2})^2 = 0.76969925830451198120784590481093\dots$:

$$
\sum_{p \equiv 2,3 \pmod{5}} \arctan\bigl(p^{-\sigma^{\ast}}\bigr) = \frac{\pi}{2} - \arctan\kappa.
$$

### 2.2 Three-Regime Stieltjes Classification & Root-Gap Obstruction (Theorem 2.4 & Proposition 2.5)

For $G(w) = \mathrm{e}^{\lambda w}\prod _ {j\ge 1}(1 + w/w _ j)$ ($\lambda \ge 0$, $w _ j > 0$, $\sum w _ j^{-1} < \infty$) with Stieltjes logarithmic derivative $S _ G(u) = \frac{\mathrm{d}}{\mathrm{d}u}\log G(u^2) = 4\lambda u + \sum _ {j\ge 1}\frac{2u}{w _ j + u^2}$:

* **Regime I** ($\lambda = 0$, $\sum w _ j^{-1/2} < \infty$, including all polynomials and upper order $\rho(G) < 1/2$): $\lim _ {u\to\infty} S _ G(u) = 0$ and $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ inside $\mathrm{Re}(s) > 1$; in particular $(1 - \beta\partial _ z^2)\Xi _ {\mathrm{DH}}(z)$ has infinitely many non-real zeros for every $\beta > 0$.
* **Regime II** ($\lambda = 0$, $\lim _ {u\to\infty} S _ G(u) = \tau \in (0,\infty){}$): $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ for all $0 < \tau < 2\sigma^{\ast} - 1$, whereas $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ for all $\tau > 2\sigma^{\ast} - 1$.
* **Regime III:** Whenever $\lambda > 0$ or lower order $\underline{\rho}(G) > 1/2$, one has $\lim _ {u\to\infty} S _ G(u) = +\infty$. Moreover, $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ whenever $\lambda > 0$, or whenever $\lambda = 0$ and $w _ {j+1} - w _ j = o(\sqrt{w _ j})$ as $j \to \infty$ (via disjoint-interval subproduct extraction and Rolle–Hurwitz monotonicity). Conversely, for every $\rho \in (1/2, 1)$, setting $R _ 1 := 16$ and $R _ k := 2^{2^{R _ {k-1}}}$ ($k \ge 2$) and placing $M _ k = \lfloor R _ k^{2\rho}\rfloor$ roots at $\sqrt{w} = R _ k$ yields upper order $\rho(G) = \rho$, lower order $\underline{\rho}(G) = 0$, $\liminf _ {u\to\infty} S _ G(u) = 0$, and infinitely many non-real zeros in $\mathrm{Re}(s) > 1$.

### 2.3 Global vs. Asymptotic Threshold Separation (Theorem 3.2)

Comparing $G _ {\mathrm{cosh},h}(w) = \cosh(2h\sqrt{w})$ and Pólya's hyperbolic sinc multiplier $G _ {\mathrm{sinc},h}(w) = \frac{\sinh(2h\sqrt{w})}{2h\sqrt{w}}$ of identical exponential type $\tau = 2h$, 384-bit FLINT/Arb Taylor-jet integration with a 512-arc Cauchy bound at $h = \frac{311}{500}$ and a 64-ball parametric Interval–Newton cover across $[\frac{310}{500}, \frac{311}{500}] \supset [\sigma^{\ast} - \frac{1}{2}, \frac{311}{500}]$ near $\rho _ {2442}$ proves unconditionally without any semigroup assumption:

$$
h_{\mathrm{global}}^{\ast}\bigl(G_{\mathrm{cosh}}; \Xi_{\mathrm{DH}}\bigr) = \sigma^{\ast} - \frac{1}{2} < \frac{311}{500} < h_{\mathrm{global}}^{\ast}\bigl(G_{\mathrm{sinc}}; \Xi_{\mathrm{DH}}\bigr).
$$

### 2.4 Strict Unconditional Bounds on $\Lambda _ {\mathrm{DH}}$ (Theorem 4.4 & Theorem 5.2)

Combining shifted Gauss–Weierstrass asymptotic simplicity ($\theta = \frac{\pi t}{8}$) and de Bruijn's equilateral-hyperbola contraction on $\lvert\mathrm{Re}(z)\rvert \le 2\gamma _ 0(t _ 1)$ with C-Arb Interval–Newton certificates at $t = \frac{1}{5}$ ($\rho _ {240}$ at 640-bit and $\rho _ {531}$ at 1350-bit) and at $t = \frac{1153}{5000} = 0.2306$ ($\rho _ {240}$ at 704-bit), together with Platt–Trudgian's bound $\Lambda _ \zeta \le \frac{1}{5}$, establishes:

$$
\Lambda_\zeta \le \frac{1}{5} < \frac{1153}{5000} < \Lambda_{\mathrm{DH}} < 2\left(\sigma^{\ast} - \frac{1}{2}\right)^2.
$$

## 3. Scope Boundaries & Anti-Overclaim Rules

* **Normalization convention:** Under $\Xi(z) = \xi(\frac{1}{2} + \frac{\mathrm{i}z}{2})$ and $\partial _ t H _ t = -\partial _ z^2 H _ t$, $\Lambda$ is $4\times$ larger than under the unscaled variable $\widetilde{\Xi}(w) = \xi(\frac{1}{2} + \mathrm{i}w)$ ($\Lambda = 4\widetilde{\Lambda}$). Always compare $\Lambda _ \zeta$ and $\Lambda _ {\mathrm{DH}}$ in the same convention.
* **Exact value of $\Lambda _ {\mathrm{DH}}$ remains open:** High-precision collision tracking yields $t _ {\mathrm{col}}(\rho _ {240}) \approx 0.2306073619$ and an exploratory scan up to $\gamma = 72{,}500$ locates no zero with larger collision time, but proving $\Lambda _ {\mathrm{DH}} = t _ {\mathrm{col}}(\rho _ {240})$ unconditionally requires a finite two-dimensional argument-principle box cover up to $\gamma _ 0(t)$.
