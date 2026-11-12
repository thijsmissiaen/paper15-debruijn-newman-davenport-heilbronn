---
name: debruijn-newman-davenport-heilbronn
description: Canonical definitions, sufficient Stieltjes Polya-Schur shift conditions, FLINT/Arb interval ball certificates, and scope boundaries for Paper XV on the de Bruijn-Newman constant and Polya-Schur multipliers for the Davenport-Heilbronn function.
---

# Skill: Pólya–Schur Shift Thresholds and Certified $\Lambda _ {\mathrm{DH}}$ Bounds (Paper XV)

## 1. Canonical Normalization & Gauss–Weierstrass Integral

For the Davenport–Heilbronn combination $\mathcal{D} _ {\mathrm{DH}}(s) = \frac{1 - \mathrm{i}\kappa}{2}L(s,\chi _ 5) + \frac{1 + \mathrm{i}\kappa}{2}L(s,\overline{\chi _ 5})$ with $\kappa = \frac{\sqrt{10 - 2\sqrt{5}} - 2}{\sqrt{5} - 1}$, the completed function under the canonical de Bruijn–Newman normalization is:

$$
\xi_{\mathrm{DH}}(s) := \left(\frac{5}{\pi}\right)^{(s+1)/2}\Gamma\left(\frac{s+1}{2}\right)\mathcal{D}_{\mathrm{DH}}(s), \qquad \Xi_{\mathrm{DH}}(z) := \xi_{\mathrm{DH}}\left(\frac{1}{2} + \frac{\mathrm{i}z}{2}\right) = 2\int_0^\infty \Phi_{\mathrm{DH}}(u)\cos(zu)\thinspace\mathrm{d}u.
$$

Under backward heat evolution $H _ t^{\mathrm{DH}}(z) = \mathrm{e}^{-t\partial _ z^2}\Xi _ {\mathrm{DH}}(z) = 2\int _ 0^\infty \mathrm{e}^{tu^2}\Phi _ {\mathrm{DH}}(u)\cos(zu)\thinspace\mathrm{d}u$, writing $z = 2\gamma - 2\mathrm{i}\delta$ ($\gamma,\delta \in \mathbb{R}$) and using $\xi _ {\mathrm{DH}}(\frac{1}{2}-v+\mathrm{i}\gamma) = \overline{\xi _ {\mathrm{DH}}(\frac{1}{2}+v+\mathrm{i}\gamma)}$ yields:

$$
H_t^{\mathrm{DH}}(2\gamma - 2\mathrm{i}\delta) = A_t(\gamma,\delta) + \overline{A_t(\gamma,-\delta)}, \qquad A_t(\gamma,\delta) := \frac{1}{\sqrt{\pi t}}\int_0^\infty \mathrm{e}^{-(v-\delta)^2/t}\thinspace\xi_{\mathrm{DH}}\left(\frac{1}{2} + v + \mathrm{i}\gamma\right)\mathrm{d}v.
$$

## 2. Certified Unconditional Results

### 2.1 Supremum Enclosure & Strict Zero Inequality (Lemma 3.1)

Solving Bombieri and Ghosh's prime-angle equation via 320-bit Arb Möbius inversion and the Intermediate Value Theorem gives $\sigma^{\ast} \in [1.12036249818332506180145892226464 \pm 10^{-30}]$, $h _ {\infty}^{\ast} = \sigma^{\ast} - \frac{1}{2} = 0.62036249818332506180145892226464\dots$, $2\sigma^{\ast} - 1 = 1.24072499636665012360291784452928\dots$, and $2(\sigma^{\ast} - \frac{1}{2})^2 = 0.76969925830451198120784590481093\dots$:

$$
\sum_{p \equiv 2,3 \pmod{5}} \arctan\bigl(p^{-\sigma^{\ast}}\bigr) = \frac{\pi}{2} - \arctan\kappa.
$$

Every individual zero $\rho = \sigma + \mathrm{i}\gamma$ of $\mathcal{D} _ {\mathrm{DH}}(s)$ satisfies $\sigma < \sigma^{\ast}$ strictly (via opposite Euler-factor signs at $\gamma = 0$ and $\frac{\log 3}{\log 2} \notin \mathbb{Q}$ for $\gamma \ne 0$).

### 2.2 Sufficient Stieltjes Shift Conditions & Lacunary Obstruction (Theorem 2.6 & Proposition 2.7)

For $G(w) = \mathrm{e}^{\lambda w}\prod _ {j\ge 1}(1 + w/w _ j)$ ($\lambda \ge 0$, $w _ j > 0$, $\sum w _ j^{-1} < \infty$) with Stieltjes logarithmic derivative $S _ G(u) = \frac{\mathrm{d}}{\mathrm{d}u}\log G(u^2) = 4\lambda u + \sum _ {j\ge 1}\frac{2u}{w _ j + u^2}$ (generalizing Dobner 2021 and using single-zero syndetic translation to avoid Righetti's 2017 non-density phenomenon):

* **Zero shift** ($\lambda = 0$, $\sum w _ j^{-1/2} < \infty$, including all polynomials and upper order $\rho(G) < 1/2$): $\lim _ {u\to\infty} S _ G(u) = 0$ and $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ inside $\lvert\mathrm{Im}(z)\rvert > 1$; in particular $(1 - \beta\partial _ z^2)\Xi _ {\mathrm{DH}}(z)$ has infinitely many non-real zeros for every $\beta > 0$.
* **Finite shift** ($\lambda = 0$, $\lim _ {u\to\infty} S _ G(u) = \tau \in (0,\infty){}$): $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; T) \gg _ G T$ for all $0 < \tau < 2\sigma^{\ast} - 1$, whereas $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ for all $\tau > 2\sigma^{\ast} - 1$.
* **Divergent shift & root-gap regularity:** Whenever $\lambda > 0$ or lower order $\underline{\rho}(G) > 1/2$, one has $\lim _ {u\to\infty} S _ G(u) = +\infty$. Moreover, $N _ {\mathrm{off}}(G(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ whenever $\lambda > 0$, or whenever $\lambda = 0$ and $w _ {j+1} - w _ j = o(\sqrt{w _ j})$ as $j \to \infty$. Conversely, for every $\rho \in (1/2, 1)$, setting $R _ 1 := 16$ and $R _ k := 2^{2^{R _ {k-1}}}$ ($k \ge 2$) and placing $M _ k = \lfloor R _ k^{2\rho}\rfloor$ roots at $\sqrt{w} = R _ k$ yields upper order $\rho(G) = \rho$, lower order $\underline{\rho}(G) = 0$, $\liminf _ {u\to\infty} S _ G(u) = 0$, and infinitely many non-real zeros in $\lvert\mathrm{Im}(z)\rvert > 1$.

### 2.3 Global Doubling Bound, Threshold Separation, and Conditional Critical Boundary (Theorem 3.2 & Proposition 3.3)

Comparing $G _ {\mathrm{cosh},h}(w) = \cosh(2h\sqrt{w})$ and Pólya's hyperbolic sinc multiplier $G _ {\mathrm{sinc},h}(w) = \frac{\sinh(2h\sqrt{w})}{2h\sqrt{w}}$ of identical exponential type $\tau = 2h$, a 384-bit FLINT/Arb 360-ball singularity-free Interval–Newton continuation across every $h \in [0, \frac{5}{8}]$ near $\rho _ {44838}$ (plus radius $10^{-30}$ point enclosures at $h = \frac{5}{8} = 0.625$ and $h = \frac{625121}{10^6} = 0.625121$) together with a computer-free interval-splitting doubling proof ($\int _ {-h}^h = \int _ {-h}^0 + \int _ 0^h$, noting that $\lbrace G _ {\mathrm{sinc},h}\rbrace _ {h>0}$ is not a semigroup in $h$) prove independently of Theorem 2.6:

$$
h_{\mathrm{global}}^{\ast}\bigl(G_{\mathrm{cosh}}; \Xi_{\mathrm{DH}}\bigr) = \sigma^{\ast} - \frac{1}{2} < \frac{5}{8} < h_{\mathrm{global}}^{\ast}\bigl(G_{\mathrm{sinc}}; \Xi_{\mathrm{DH}}\bigr) \le 2\sigma^{\ast} - 1 = 1.240724996\dots
$$

At the critical asymptotic boundary $h _ {\ast} = \sigma^{\ast} - \frac{1}{2}$ ($\tau = 2\sigma^{\ast} - 1$), ball Interval–Newton enclosures on $H _ {\ast} = [\sigma^{\ast} - \frac{1}{2}]$ near $\rho _ {2442}$ and $\rho _ {44838}$ certify unconditionally $N _ {\mathrm{off}}(G _ {\mathrm{sinc},h _ {\ast}}(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) \ge 8$, while conditionally under RH for $\zeta(s)$ and GRH for $L(s, \chi _ {5,\mathrm{quad}})$, Proposition 3.3 proves $8 \le N _ {\mathrm{off}}(G _ {\mathrm{sinc},h _ {\ast}}(-\partial _ z^2)\Xi _ {\mathrm{DH}}; \infty) < \infty$ because the Borel–Laplace secondary tail $\mathcal{O}((\log\gamma)^{-1})$ is strictly dominated by the conditional prime-angle lower bound $\lvert\mathcal{D} _ {\mathrm{DH}}(s(z)+h _ {\ast})\rvert \gg (\log\gamma)^{-1/3}$. Unconditionally, the Vinogradov–Korobov zero-free region and Baker–Wüstholz linear forms in logarithms leave an exponential gap against $\mathcal{O}((\log\gamma)^{-1})$ (Remark 3.4).

### 2.4 Unconditional Bounds on $\Lambda _ {\mathrm{DH}}$ & Dobner RS $\ell^1$ Barrier (Corollary 4.4, Remark 4.5, & Theorem 5.2)

Adapting the large-height simplicity and finite-box hyperbola argument of Ki, Kim, and Lee (2009) yields $\Lambda _ {\mathrm{DH}} < 2(\sigma^{\ast} - \frac{1}{2})^2$ (whereas Dobner's Riemann–Siegel-truncated modified Dirichlet polynomial of length $N(\gamma) = \lfloor\sqrt{5\gamma/(2\pi)}\rfloor$ suffers heat amplification $\mathrm{e}^{\frac{t}{4}\log^2 n}$, inflating the $\ell^1$ dominance height from $\gamma _ {\ell^1}(0.50) \approx 1.52 \times 10^4$ to $\gamma _ {\ell^1}(0.25) \approx 4.18 \times 10^8$; Remark 4.5). Meanwhile, a 512-bit FLINT/Arb Interval–Newton enclosure (radius $10^{-32}$) at exact rational $t = \frac{28825920237}{125000000000} = 0.230607361896$ near $\rho _ {240}$ establishes:

$$
\frac{28825920237}{125000000000} < \Lambda_{\mathrm{DH}} < 2\left(\sigma^{\ast} - \frac{1}{2}\right)^2.
$$

## 3. Scope Boundaries & Anti-Overclaim Rules

* **Sufficient conditions, not an exhaustive partition:** Oscillating $S _ G(u)$ crossing $2\sigma^{\ast} - 1$ without root-gap regularity, the unconditional status of the critical boundary $\lim _ {u\to\infty} S _ G(u) = 2\sigma^{\ast} - 1$ (finiteness for $G _ {\mathrm{sinc},\sigma^{\ast}-1/2}$ is conditional on RH for $\zeta(s)$ and GRH for $L(s, \chi _ {5,\mathrm{quad}})$ in Proposition 3.3), and clustered multipliers with $\underline{\rho}(G) > 1/2$ violating $w _ {j+1} - w _ j = o(\sqrt{w _ j})$ remain open (Remarks 2.8 and 3.4).
* **Normalization & open exact value:** Under $\Xi(z) = \xi(\frac{1}{2} + \frac{\mathrm{i}z}{2})$ and $\partial _ t H _ t = -\partial _ z^2 H _ t$, $\Lambda = 4\widetilde{\Lambda}$ (so Platt–Trudgian's $\Lambda _ \zeta \le \frac{1}{5}$ compares directly with $\Lambda _ {\mathrm{DH}} > \frac{28825920237}{125000000000}$). High-precision scalar Newton tracking gives $t _ {\mathrm{col}}(\rho _ {240}) \approx 0.2306073618964546$ (matching all 12 decimal digits of $\frac{28825920237}{125000000000}$) and the curve $4(\sigma^{\ast} - \frac{1}{2})/\log(5\gamma/2\pi)$ is only an unproved leading-order heuristic; proving $\Lambda _ {\mathrm{DH}} = t _ {\mathrm{col}}(\rho _ {240})$ remains open.
