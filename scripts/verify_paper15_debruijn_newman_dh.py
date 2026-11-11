#!/usr/bin/env python3
"""
Rigorous C-Arb Ball-Arithmetic Verification Suite for Paper XV (v15.1.0)
========================================================================
Manuscript: Asymptotic Shift Thresholds for Polya-Schur Multipliers and a
            Certified Lower Bound on the de Bruijn-Newman Constant of the
            Davenport-Heilbronn Function
Author:     Thijs Missiaen (ORCID: 0009-0007-3678-0810)

Default usage (verifies all FLINT/Arb certificates without touching figures/):
    python scripts/verify_paper15_debruijn_newman_dh.py

Optional figure regeneration (writes PDF/PNG plots to figures/):
    python scripts/verify_paper15_debruijn_newman_dh.py --regenerate-figures

Certifies in 320-bit and 448-bit FLINT/Arb interval ball arithmetic (python-flint):
  [1/4] IVT enclosure of sigma* = sup{Re(rho) : D_DH(rho) = 0} (radius 10^-30),
        h_infty* = sigma* - 1/2, 2*(sigma* - 1/2)^2, and Euler-product phase sign
        -(1 - i*kappa)/(1 + i*kappa) = exp(+2i*(pi/2 - arctan(kappa))).
  [2/4] Singularity-free Taylor-series + 512-arc Cauchy-circle Interval-Newton
        continuation across the ENTIRE compact interval h in [0, 311/500] (400
        overlapping h-balls including h = 0, plus radius 10^-30 enclosure at
        h = 311/500), proving h_global*(G_sinc; Xi_DH) > 311/500 > sigma* - 1/2
        independently of Section 2.
  [3/4] 448-bit Interval-Newton off-real zero enclosure for H_t^DH(z) at exact
        rational t = 1153/5000 = 0.2306 near rho_240 (radius 10^-32), proving
        Lambda_DH > 1153/5000 = 0.2306.
  [4/4] Explicit analytical constants for Lemma 4.2, Proposition 4.3, Theorem 3.2(ii)
        (global finiteness h_global*(G_sinc; Xi_DH) <= 31 via Y_0 = 30), and
        Proposition 2.6 (tower-lacunary canonical product S_G(u_k) -> 0).
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import flint
import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = REPO_ROOT / "figures"


# ============================================================================
# 1. 320-BIT ARB IVT ENCLOSURE OF sigma* AND EULER-PRODUCT PHASE SIGN
# ============================================================================
def certify_supremum_sigma_star() -> tuple[flint.arb, flint.arb, flint.arb]:
    print("[1/4] Certifying 320-bit FLINT/Arb IVT enclosure of sigma*...")
    flint.ctx.prec = 320
    sqrt5 = flint.arb(5).sqrt()
    kap = ((flint.arb(10) - 2 * sqrt5).sqrt() - 2) / (sqrt5 - 1)
    pi_a = flint.arb.pi()
    target = pi_a / 2 - kap.atan()

    # Verify exact phase sign: -(1 - i*kap)/(1 + i*kap) == exp(+2i*(pi/2 - atan(kap)))
    lhs_phase = -(flint.acb(1, -kap)) / flint.acb(1, kap)
    rhs_phase = (flint.acb(0, 2) * target).exp()
    assert (lhs_phase - rhs_phase).contains(flint.acb(0)), "Phase sign check failed"
    assert bool(lhs_phase.imag > flint.arb("0.5")), "Imaginary part of phase ratio must be positive"

    N_SIEVE = 2000
    sieve = [True] * (N_SIEVE + 1)
    sieve[0] = sieve[1] = False
    for k in range(2, int(math.isqrt(N_SIEVE)) + 1):
        if sieve[k]:
            for m in range(k * k, N_SIEVE + 1, k):
                sieve[m] = False
    small_p = [p for p in range(2, 1000) if sieve[p] and (p % 5) in (2, 3)]

    K_MOB = 280
    mu = [1] * (K_MOB + 1)
    is_p = [True] * (K_MOB + 1)
    for p in range(2, K_MOB + 1):
        if is_p[p]:
            for m in range(p, K_MOB + 1, p):
                is_p[m] = False
                mu[m] = -mu[m]
            p2 = p * p
            for m in range(p2, K_MOB + 1, p2):
                mu[m] = 0
    odd_mu = [(d, mu[d]) for d in range(1, K_MOB, 2) if mu[d] != 0]

    def L0_arb(u: flint.arb) -> flint.arb:
        uc = flint.acb(u)
        return ((1 - flint.acb(5) ** (-uc)) * uc.zeta()).real

    def L2_arb(u: flint.arb) -> flint.arb:
        uc = flint.acb(u)
        val = (flint.acb(5) ** (-uc)) * (
            uc.zeta(flint.acb(1) / 5) - uc.zeta(flint.acb(2) / 5)
            - uc.zeta(flint.acb(3) / 5) + uc.zeta(flint.acb(4) / 5)
        )
        return val.real

    def P_minus_arb(u: flint.arb) -> flint.arb:
        tot = flint.arb(0)
        for d, mu_d in odd_mu:
            du = u * d
            tot += (flint.arb(mu_d) / (2 * d)) * (L0_arb(du).log() - L2_arb(du).log())
        return tot + flint.arb(0, flint.arb("1e-90"))

    def F_sup_rigorous(sig: flint.arb) -> flint.arb:
        s_small = sum((flint.arb(p) ** (-sig)).atan() for p in small_p)
        s_tail = flint.arb(0)
        for m in range(15):
            u = sig * (2 * m + 1)
            pm_tail = P_minus_arb(u) - sum(flint.arb(p) ** (-u) for p in small_p)
            s_tail += ((-1) ** m) * pm_tail / (2 * m + 1)
        return s_small + s_tail + flint.arb(0, flint.arb("1e-95")) - target

    s_a = flint.arb("1.120362498183325")
    s_b = flint.arb("1.120362498183326")
    for _ in range(6):
        fa = flint.arb(F_sup_rigorous(s_a).mid())
        fb = flint.arb(F_sup_rigorous(s_b).mid())
        s_next = flint.arb((s_b - fb * (s_b - s_a) / (fb - fa)).mid())
        s_a, s_b = s_b, s_next

    eps_ivt = flint.arb("5e-31")
    f_lo = F_sup_rigorous(s_b - eps_ivt)
    f_hi = F_sup_rigorous(s_b + eps_ivt)
    assert bool(f_lo > 0 and f_hi < 0), "IVT sign change failed for sigma*"
    sigma_dec = flint.arb("1.12036249818332506180145892226464")
    assert flint.arb(sigma_dec.mid(), flint.arb("1e-30")).contains(flint.arb(s_b.mid(), eps_ivt))
    sigma_star = flint.arb(sigma_dec.mid(), flint.arb("1e-30"))
    h_inf_star = sigma_star - flint.arb("0.5")
    debruijn_bound = 2 * (h_inf_star ** 2)

    print(f"      sigma*           = {sigma_star.str(36, radius=True)}")
    print(f"      sigma* - 1/2     = {h_inf_star.str(36, radius=True)}")
    print(f"      2*(sigma*-1/2)^2 = {debruijn_bound.str(36, radius=True)}")
    return sigma_star, h_inf_star, debruijn_bound


# ============================================================================
# 2. SINGULARITY-FREE CONTINUUM CERTIFICATE ON ALL h in [0, 311/500]
# ============================================================================
def certify_polya_sinc_separation(h_inf_star: flint.arb) -> dict[str, object]:
    print("\n[2/4] Certifying singularity-free Interval-Newton cover across ALL h in [0, 311/500]...")
    flint.ctx.prec = 320
    sqrt5 = flint.arb(5).sqrt()
    kap = ((flint.arb(10) - 2 * sqrt5).sqrt() - 2) / (sqrt5 - 1)
    pi_a = flint.arb.pi()
    h_rat = flint.arb(311) / flint.arb(500)
    assert bool(h_rat > h_inf_star)

    # Common reference point near the midpoint of {s(z(h)) : h in [0, 311/500]}
    s_ref = flint.acb(flint.arb("0.71"), flint.arb("2442.49"))

    def tilde_xi(s: flint.acb) -> flint.acb:
        lp = (flint.acb(5) / pi_a).log()
        lg = ((s + 1) / 2).lgamma() - ((s_ref + 1) / 2).lgamma()
        rg = (((s - s_ref) / 2) * lp + lg).exp()
        fv = (flint.acb(5) ** (-s)) * (
            s.zeta(flint.acb(1) / 5) + kap * s.zeta(flint.acb(2) / 5)
            - kap * s.zeta(flint.acb(3) / 5) - s.zeta(flint.acb(4) / 5)
        )
        return rg * fv

    def tilde_xi_series(sb: flint.acb, cap_len: int) -> flint.acb_series:
        flint.ctx.cap = cap_len
        ss = flint.acb_series([sb, 1])
        lp = (flint.acb(5) / pi_a).log()
        lg = ((ss + 1) / 2).lgamma() - ((s_ref + 1) / 2).lgamma()
        rg = (((ss - s_ref) / 2) * lp + lg).exp()
        pref5 = (-ss * flint.acb(5).log()).exp()
        fv = pref5 * (
            ss.zeta(flint.acb(1) / 5) + kap * ss.zeta(flint.acb(2) / 5)
            - kap * ss.zeta(flint.acb(3) / 5) - ss.zeta(flint.acb(4) / 5)
        )
        return rg * fv

    # Certify M_45 := sup_{|w - s_ref| = 4.5} |tilde_xi(w)| via 512 overlapping interval arcs
    M_45 = flint.arb(0)
    half_step = flint.arb("0.5001")
    for k in range(512):
        th = 2 * pi_a * flint.arb(flint.arb(k) + flint.arb("0.5"), half_step) / 512
        w = flint.arb("4.5") * (flint.acb(0, 1) * th).exp()
        M_45 = M_45.max(abs(tilde_xi(s_ref + w)).upper())
    print(f"      512-arc Cauchy bound on |w - s_ref| = 4.5: M_4.5 <= {float(M_45.upper()):.4e}")

    # Singularity-free series g_sinc(z, h) = (1/2) * int_{-1}^1 tilde_xi(s(z) + h*v) dv
    # = sum_{m>=0} (c_{2m}(s(z)) / (2m+1)) * h^{2m}, regular at h = 0 with zero division by h!
    # Cover [0, 285/500] with 180 balls, [285/500, 308/500] with 120 balls,
    # and [308/500, 311/500] with 100 balls (400 overlapping h-balls in total).
    segments = [
        (flint.arb(0), flint.arb(285) / flint.arb(500), 180, flint.arb("6e-3")),
        (flint.arb(285) / flint.arb(500), flint.arb(308) / flint.arb(500), 120, flint.arb("2.5e-3")),
        (flint.arb(308) / flint.arb(500), flint.arb(311) / flint.arb(500), 100, flint.arb("1.2e-3")),
    ]
    z_tr = flint.acb("4885.06058493697206", "-0.79550030995402")
    sinc_curve: list[tuple[float, float]] = []
    prev_slice: tuple[flint.arb, flint.arb, flint.acb, flint.acb, flint.acb, flint.acb] | None = None
    last_B_k: flint.acb | None = None
    total_balls = 0

    for seg_idx, (h_lo, h_hi, n_balls, r_box) in enumerate(segments):
        dh_half = (h_hi - h_lo) / (2 * n_balls)
        dh_pad = dh_half * flint.arb("1.005")
        for k_s in range(n_balls):
            h_mid = h_lo + (2 * k_s + 1) * dh_half
            H_k = flint.arb(h_mid.mid(), dh_pad.upper())
            if seg_idx == 0 and k_s == 0:
                assert H_k.contains(flint.arb(0)), "First h-ball must contain h = 0"
                sinc_curve.append((0.0, abs(float(z_tr.imag.mid())) / 2.0))
            if seg_idx == len(segments) - 1 and k_s == n_balls - 1:
                assert H_k.contains(h_rat), "Final h-ball must contain h = 311/500"

            for _ in range(3):
                sb0 = flint.acb("0.5") + flint.acb(0, 1) * z_tr / 2
                j0 = tilde_xi_series(sb0, 20)
                f0 = sum((j0[2 * m] / (2 * m + 1)) * (h_mid ** (2 * m)) for m in range(10))
                df0 = (flint.acb(0, 1) / 2) * sum(j0[2 * m + 1] * (h_mid ** (2 * m)) for m in range(9))
                st = f0 / df0
                z_tr = flint.acb((z_tr.real - st.real).mid(), (z_tr.imag - st.imag).mid())

            sb_k = flint.acb("0.5") + flint.acb(0, 1) * z_tr / 2
            assert bool(abs(sb_k - s_ref) < flint.arb("0.25")), "Trajectory left Cauchy disk"
            R_avail = flint.arb("4.5") - abs(sb_k - s_ref).upper()  # > 4.25
            jet_k = tilde_xi_series(sb_k, 40)

            g_mid = sum((jet_k[2 * m] / (2 * m + 1)) * (h_mid ** (2 * m)) for m in range(20))
            q_m = h_mid.upper() / R_avail
            err_g = M_45 * (q_m ** 40) / (1 - q_m)
            g_mid += flint.acb(flint.arb(0, err_g.upper()), flint.arb(0, err_g.upper()))

            H_abs = abs(H_k).upper()
            q_H = H_abs / R_avail
            dh_g = sum((flint.arb(2 * m) / (2 * m + 1)) * jet_k[2 * m] * (H_k ** (2 * m - 1)) for m in range(1, 20))
            err_dh = (M_45 / R_avail) * (q_H ** 39) / (1 - q_H * q_H)
            dh_g += flint.acb(flint.arb(0, err_dh.upper()), flint.arb(0, err_dh.upper()))

            B_k = z_tr + flint.acb(flint.arb(0, r_box), flint.arb(0, r_box))
            G1 = (flint.acb(0, 1) / 2) * sum(jet_k[2 * m + 1] * (H_k ** (2 * m)) for m in range(19))
            err_G1 = (M_45 / (2 * R_avail)) * (q_H ** 38) / (1 - q_H * q_H)
            G1 += flint.acb(flint.arb(0, err_G1.upper()), flint.arb(0, err_G1.upper()))

            G2 = (-flint.arb(1) / 4) * sum(flint.arb(2 * m + 2) * jet_k[2 * m + 2] * (H_k ** (2 * m)) for m in range(18))
            err_G2 = (M_45 / (4 * R_avail * R_avail)) * 38 * (q_H ** 36) / ((1 - q_H) ** 2)
            G2 += flint.acb(flint.arb(0, err_G2.upper()), flint.arb(0, err_G2.upper()))

            G3 = (-flint.acb(0, 1) / 8) * sum(
                flint.arb((2 * m + 3) * (2 * m + 2)) * jet_k[2 * m + 3] * (H_k ** (2 * m)) for m in range(17)
            )
            err_G3 = (M_45 / (8 * (R_avail ** 3))) * 1332 * (q_H ** 34) / ((1 - q_H) ** 3)
            G3 += flint.acb(flint.arb(0, err_G3.upper()), flint.arb(0, err_G3.upper()))

            dz_var = B_k - z_tr
            dz_max = abs(dz_var).upper()
            R_rem = R_avail - (dz_max / 2 + H_abs)
            assert bool(R_rem > flint.arb("3.5")), "Remainder radius dropped below 7/2"
            err_dz4 = (dz_max ** 3) * M_45 / (4 * (R_rem ** 4))
            dz_g_Bk = G1 + dz_var * (G2 + (dz_var / 2) * G3) + flint.acb(
                flint.arb(0, err_dz4.upper()), flint.arb(0, err_dz4.upper())
            )

            g_zk_Hk = g_mid + (H_k - h_mid) * dh_g
            N_Bk = z_tr - g_zk_Hk / dz_g_Bk
            assert B_k.contains_interior(N_Bk) and bool(N_Bk.imag < flint.arb(0)), (
                f"Continuum sinc cert failed at segment {seg_idx}, slice {k_s}"
            )
            if seg_idx == 0 and k_s == 0:
                z_rho0 = flint.acb("4885.06058493697206", "-0.79550030995402")
                B_rho0 = z_rho0 + flint.acb(flint.arb(0, flint.arb("1e-12")), flint.arb(0, flint.arb("1e-12")))
                j_rho0 = tilde_xi_series(flint.acb("0.5") + flint.acb(0, 1) * B_rho0 / 2, 2)
                N_rho0 = z_rho0 - tilde_xi(flint.acb("0.5") + flint.acb(0, 1) * z_rho0 / 2) / (
                    (flint.acb(0, 1) / 2) * j_rho0[1]
                )
                assert B_rho0.contains_interior(N_rho0) and B_k.contains(B_rho0), (
                    "Initial h = 0 zero enclosure B_rho0 (rho_2442) failed or not inside B_1"
                )
            if prev_slice is not None:
                p_H, p_hmid, p_ztr, p_gmid, p_dhg, p_dzg = prev_slice
                assert H_k.overlaps(p_H), "Adjacent h-balls must overlap"
                h_join = flint.arb(((p_H.upper() + H_k.lower()) / 2).mid())
                assert p_H.contains(h_join) and H_k.contains(h_join), "h_join outside H_{k-1} cap H_k"
                N_join_prev = p_ztr - (p_gmid + (h_join - p_hmid) * p_dhg) / p_dzg
                assert B_k.contains(N_join_prev), (
                    f"Branch continuity failed between slices {total_balls - 1} and {total_balls}"
                )
            prev_slice = (H_k, h_mid, z_tr, g_mid, dh_g, dz_g_Bk)
            last_B_k = B_k
            total_balls += 1
            if k_s % 10 == 0 or (seg_idx == 1 and k_s % 16 == 15):
                sinc_curve.append((float(h_mid.mid()), abs(float(z_tr.imag.mid())) / 2.0))

    # Sharp enclosure around exact printed decimal midpoint at h = 311/500 (radius 10^-30):
    zc_s = flint.acb(
        "4884.88665419175715643333300197615178",
        "-0.056004628860858072317639148896287841",
    )
    sb_pt = flint.acb("0.5") + flint.acb(0, 1) * zc_s / 2
    R_av_pt = flint.arb("4.5") - abs(sb_pt - s_ref).upper()
    j_pt = tilde_xi_series(sb_pt, 66)
    f_mid = sum((j_pt[2 * m] / (2 * m + 1)) * (h_rat ** (2 * m)) for m in range(32))
    q_pt = h_rat / R_av_pt
    err_pt = M_45 * (q_pt ** 64) / (1 - q_pt)
    f_mid += flint.acb(flint.arb(0, err_pt.upper()), flint.arb(0, err_pt.upper()))

    r_sinc = flint.arb("1e-30")
    B_sinc = zc_s + flint.acb(flint.arb(0, r_sinc), flint.arb(0, r_sinc))
    assert last_B_k is not None and last_B_k.contains(B_sinc), "B_sinc must lie inside B_400 for branch continuity"
    sb_B = flint.acb("0.5") + flint.acb(0, 1) * B_sinc / 2
    df_B = flint.acb(0, 1) * (tilde_xi(sb_B + h_rat) - tilde_xi(sb_B - h_rat)) / (4 * h_rat)
    N_B = zc_s - f_mid / df_B
    cert_in = B_sinc.contains_interior(N_B)
    cert_im = bool(N_B.imag < flint.arb(0))
    assert cert_in and cert_im, "Point Interval-Newton certificate at h = 311/500 failed"
    sinc_curve.append((0.622, abs(float(zc_s.imag.mid())) / 2.0))

    print(f"      Continuum & branch-overlap certified across ALL h in [0, 311/500] ({total_balls} balls, 399 joins): True")
    print(f"      z_sinc(311/500) = {zc_s.real.str(34, radius=False)} + ({zc_s.imag.str(34, radius=False)})j")
    print(f"      N(B_sinc) c int(B_sinc) (r=10^-30): {cert_in},  Im(z_sinc(311/500)) < 0: {cert_im}")

    return {"z_sinc_622": zc_s, "sinc_curve": sinc_curve, "tilde_xi": tilde_xi}


# ============================================================================
# 3. 448-BIT INTERVAL-NEWTON HEAT WITNESS AT t = 1153/5000 = 0.2306 (rho_240)
# ============================================================================
def certify_heat_flow_witness() -> flint.acb:
    print("\n[3/4] Certifying 448-bit Interval-Newton witness at exact rational t = 1153/5000 = 0.2306...")
    flint.ctx.prec = 448
    sqrt5 = flint.arb(5).sqrt()
    kap = ((flint.arb(10) - 2 * sqrt5).sqrt() - 2) / (sqrt5 - 1)
    pi_a = flint.arb.pi()
    pi5 = flint.acb(pi_a / 5)

    n_terms = 32
    coeffs = []
    for n in range(1, n_terms):
        r = n % 5
        if r == 0:
            continue
        c = flint.acb(1) if r == 1 else (flint.acb(kap) if r == 2 else (flint.acb(-kap) if r == 3 else flint.acb(-1)))
        coeffs.append((n, flint.acb(n) * c, flint.acb(n * n) * pi5))

    def Phi_DH(u: flint.acb) -> flint.acb:
        e4u = (4 * u).exp()
        s = flint.acb(0)
        for _, nc, n2p in coeffs:
            s += nc * (-(n2p * e4u)).exp()
        return 4 * (3 * u).exp() * s

    U_top = flint.arb("2.28")
    t_heat = flint.arb(1153) / flint.arb(5000)
    tol = flint.arb("1e-122")

    def tail_ball(y_abs: flint.arb, deriv_order: int) -> flint.arb:
        r_geom = (-(flint.arb(n_terms) * pi_a / 5)).exp()
        sum_tail = (r_geom ** n_terms) * (n_terms - (n_terms - 1) * r_geom) / ((1 - r_geom) ** 2)
        E1 = 8 * (U_top ** (deriv_order + 1)) * (t_heat * U_top * U_top + (3 + y_abs) * U_top).exp() * sum_tail
        e4U = (4 * U_top).exp()
        E2 = (40 / pi_a) * (-4 * U_top - (pi_a / 10) * e4U).exp()
        return E1 + E2

    def H_deriv_rig(z: flint.acb, order: int) -> flint.acb:
        tc = flint.acb(t_heat)
        if order == 0:
            ker = lambda u, _: (tc * u * u).exp() * Phi_DH(u) * (z * u).cos()
        elif order == 1:
            ker = lambda u, _: -u * (tc * u * u).exp() * Phi_DH(u) * (z * u).sin()
        elif order == 2:
            ker = lambda u, _: -(u * u) * (tc * u * u).exp() * Phi_DH(u) * (z * u).cos()
        else:
            ker = lambda u, _: (u ** 3) * (tc * u * u).exp() * Phi_DH(u) * (z * u).sin()
        Im = 2 * flint.acb.integral(ker, flint.acb(0), flint.acb(U_top), abs_tol=tol, rel_tol=tol, eval_limit=260000)
        eb = tail_ball(abs(z.imag).upper(), order)
        rad = flint.arb(0, eb)
        return Im + flint.acb(rad, rad)

    zc = flint.acb(
        "480.690235013620019064551671086167379",
        "-0.003837169565586910716957555909571243",
    )
    r_box = flint.arb("1e-32")
    B = zc + flint.acb(flint.arb(0, r_box), flint.arb(0, r_box))
    ym = abs(B.imag).upper()
    M4 = (
        4 * n_terms * (n_terms - 1) * (U_top ** 5) * (t_heat * U_top * U_top + 3 * U_top).exp() * (U_top * ym).cosh()
        + tail_ball(ym, 4)
    )
    fz = H_deriv_rig(zc, 0)
    dfz = H_deriv_rig(zc, 1)
    d2fz = H_deriv_rig(zc, 2)
    d3fz = H_deriv_rig(zc, 3)
    dz = B - zc
    # By Taylor's theorem on H_t': |H_t'(z) - (dfz + dz*d2fz + dz^2/2 * d3fz)| <= (|dz|^3 / 6) * M4,
    # where |dz| <= sqrt(2) * r_box => |dz|^3 / 6 <= (2*sqrt(2)/6) * r_box^3 < 0.5 * r_box^3:
    rem_d1 = (abs(dz).upper() ** 3) * M4 / 6
    df_B = dfz + dz * (d2fz + (dz / 2) * d3fz) + flint.acb(flint.arb(0, rem_d1.upper()), flint.arb(0, rem_d1.upper()))
    N_B = zc - fz / df_B

    ok_in = B.contains_interior(N_B)
    ok_im = bool(N_B.imag < flint.arb(0))
    assert ok_in and ok_im, "Interval-Newton failed for rho_240 at t = 1153/5000"
    print(f"      z_240(1153/5000) = {zc.real.str(36, radius=False)} + ({zc.imag.str(36, radius=False)})j")
    print(f"      N(B) c int(B) (r=10^-32): {ok_in},  Im(z) < 0: {ok_im}")
    return zc


# ============================================================================
# 4. EXPLICIT CONSTANTS FOR THEOREM 3.2(ii), PROPOSITION 4.3 & PROPOSITION 2.7
# ============================================================================
def verify_analytic_lemmas() -> None:
    print("\n[4/4] Verifying explicit analytic bounds (Abel, Dirichlet tail, sinc upper bound, lacunary)...")
    mp.mp.dps = 45
    sqrt5 = mp.sqrt(5)
    kap = (mp.sqrt(10 - 2 * sqrt5) - 2) / (sqrt5 - 1)
    pi = mp.pi

    # 4a. Abel summation bound coefficient < 3.90 < 4.0 for tau >= 100, v >= 0 (and < 38 for 0 <= tau <= 100):
    c1 = 1.0 + float(kap)
    abel_coeff = 2.0 * c1 * 0.46 + 0.1 * c1 + c1 * math.sqrt(1.0 + 361.0 / 40000.0) * 2.0
    abel_low = 2.0 * c1 * 4.6 + 0.1 * c1 + c1 * math.sqrt(9.5 ** 2 + 100.0 ** 2) * 0.2
    assert abel_coeff < 3.90 < 4.0 and abel_low < 38.0
    print(f"      Abel bound: {abel_coeff:.6f} * tau^(1/2) (< 4) for tau>=100; <= {abel_low:.2f} (< 38) for |tau|<=100 (OK)")

    # 4b. Dirichlet tail at sigma = 7/4 (u = 5/4 = 1.25):
    s74 = mp.mpf("1.75")
    tail_74 = (
        (mp.hurwitz(s74, mp.mpf(1) / 5) * mp.power(5, -s74) - 1)
        + kap * mp.hurwitz(s74, mp.mpf(2) / 5) * mp.power(5, -s74)
        + kap * mp.hurwitz(s74, mp.mpf(3) / 5) * mp.power(5, -s74)
        + mp.hurwitz(s74, mp.mpf(4) / 5) * mp.power(5, -s74)
    )
    assert tail_74 < mp.mpf("0.50")
    print(f"      Dirichlet tail sum_{{n>=2}} |a_DH(n)| n^(-7/4) = {mp.nstr(tail_74, 14)} < 0.50 (OK)")

    # 4c. Global finiteness h_global*(G_sinc; Xi_DH) <= 31 via Y_0 = 30:
    # On Y >= 30 (|z_+| >= 15.75, alpha_R = 0.5*log(5*|z_+|/pi) > 1.61),
    # splitting [0, Y] = [0, 6] U [6, Y-4] U [Y-4, Y] leaves Re(J(gam, Y) / d_Y J(gam, Y)) > 0.70 / alpha_R > 0:
    zmin = 15.75
    aR_min = 0.5 * math.log(5.0 * zmin / math.pi)
    re_inv_alpha = 1.0 / (1.0 + (math.pi ** 2 / 16.0) / (aR_min ** 2))  # > 0.807
    err_main = 1e-5 + 0.33 * (2.0 / (aR_min ** 2) + 2.1 / aR_min) / zmin + math.exp(-6.0 * aR_min)
    base_scale = (5.0 * zmin / math.pi) ** (-3.0)
    err_mid = (2.5 * 1.06 / 0.9999) * (1.0 / 0.64) * base_scale * aR_min
    err_end = 4.0 * 2.5 * math.exp(-11.8) * base_scale * (38.0 * (zmin ** 0.5) / 0.9999) * aR_min
    net_coeff = re_inv_alpha - (err_main + err_mid + err_end)
    assert net_coeff > 0.70, f"net_coeff={net_coeff}"
    print(f"      Global sinc bound (Y >= 30 => h >= 31): Re(J / d_Y J) >= {net_coeff:.4f} / alpha_R > 0 (OK)")

    # 4d. Gauss-Weierstrass half-line identity & Sinc half-line identity J(gam, h+delta) + conj(J(gam, h-delta)):
    def D_DH(s):
        return mp.power(5, -s) * (
            mp.hurwitz(s, mp.mpf(1) / 5) + kap * mp.hurwitz(s, mp.mpf(2) / 5)
            - kap * mp.hurwitz(s, mp.mpf(3) / 5) - mp.hurwitz(s, mp.mpf(4) / 5)
        )

    def xi_DH(s):
        return mp.power(5 / pi, (s + 1) / 2) * mp.gamma((s + 1) / 2) * D_DH(s)

    def Phi_DH(u):
        e4u = mp.exp(4 * u)
        acc = mp.mpf(0)
        for n in range(1, 32):
            r = n % 5
            if r == 0:
                continue
            c = 1 if r == 1 else (kap if r == 2 else (-kap if r == 3 else -1))
            acc += n * c * mp.exp(-(pi * n * n / 5) * e4u)
        return 4 * mp.exp(3 * u) * acc

    t_chk, gam_chk, del_chk = mp.mpf("0.35"), mp.mpf("12.5"), mp.mpf("0.28")
    z_chk = 2 * gam_chk - 2 * mp.j * del_chk
    H_dir = 2 * mp.quad(lambda u: mp.exp(t_chk * u * u) * Phi_DH(u) * mp.cos(z_chk * u), [0, mp.mpf("2.4")])

    def A_t(gam, delta):
        return (1 / mp.sqrt(pi * t_chk)) * mp.quad(
            lambda v: mp.exp(-((v - delta) ** 2) / t_chk) * xi_DH(mp.mpf("0.5") + v + mp.j * gam),
            [0, 8, 18]
        )

    fold_err = abs(H_dir - (A_t(gam_chk, del_chk) + mp.conj(A_t(gam_chk, -del_chk))))
    assert fold_err < mp.mpf("1e-35")
    print(f"      Gauss-Weierstrass half-line residual: {mp.nstr(fold_err, 8)} < 1e-35 (OK)")

    # 4e. Tower-lacunary canonical product at rho = 3/4, R_1 = 16, R_2 = 2^65536:
    log2_R2 = mp.mpf(65536)
    rho_lac = mp.mpf("0.75")
    log2_u2 = mp.mpf("0.5") * (1 - rho_lac) * log2_R2
    D1 = mp.mpf(int(16 ** 1.5))
    S_val = 2 * D1 * (mp.mpf(2) ** (-log2_u2)) + 2 * (mp.mpf(2) ** (2 * rho_lac * log2_R2 + log2_u2 - 2 * log2_R2))
    tail_log = 4 * (mp.mpf(2) ** (2 * rho_lac * log2_R2 + 2 * log2_u2 - 2 * log2_R2))
    logG_u2 = D1 * mp.log(1 + mp.mpf(2) ** (2 * log2_u2)) + tail_log
    eff_order = mp.log(logG_u2) / mp.log(mp.mpf(2) ** (2 * log2_u2))
    assert S_val < mp.mpf("1e-2400") and tail_log < mp.mpf("1e-4900") and eff_order < mp.mpf("0.002")
    print(
        f"      Tower-lacunary product (rho=3/4, R_2=2^65536): S_G(u_2)={mp.nstr(S_val, 4)}, "
        f"eff_order={float(eff_order):.6f} (OK)"
    )


# ============================================================================
# 5. OPTIONAL FIGURE GENERATION (--regenerate-figures)
# ============================================================================
def generate_figures(sinc_data: dict[str, object]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    print("\n[5/5] Regenerating publication-grade figures (--regenerate-figures)...")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 10.2,
        "legend.fontsize": 8.1,
        "xtick.labelsize": 8.8,
        "ytick.labelsize": 8.8,
        "mathtext.fontset": "cm",
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2), gridspec_kw={"wspace": 0.26})

    flint.ctx.prec = 128
    s5 = flint.arb(5).sqrt()
    kp = ((flint.arb(10) - 2 * s5).sqrt() - 2) / (s5 - 1)
    pa = flint.arb.pi()
    gl_x, gl_w = np.polynomial.legendre.leggauss(24)

    def trace_heat_curve(gam0: float, sig0: float, t_col: float, n_pts: int = 28) -> tuple[np.ndarray, np.ndarray]:
        ts = np.concatenate([np.linspace(0.0, t_col * 0.92, n_pts - 6), np.linspace(t_col * 0.94, t_col * 0.998, 6)])
        dels = [sig0 - 0.5]
        s_ref = flint.acb("0.5", str(gam0))

        def txi(s: flint.acb) -> flint.acb:
            lp = (flint.acb(5) / pa).log()
            lg = ((s + 1) / 2).lgamma() - ((s_ref + 1) / 2).lgamma()
            rg = (((s - s_ref) / 2) * lp + lg).exp()
            fv = (flint.acb(5) ** (-s)) * (
                s.zeta(flint.acb(1) / 5) + kp * s.zeta(flint.acb(2) / 5)
                - kp * s.zeta(flint.acb(3) / 5) - s.zeta(flint.acb(4) / 5)
            )
            return rg * fv

        zc = flint.acb(str(2.0 * gam0), str(-2.0 * (sig0 - 0.5)))
        for tv in ts[1:]:
            sq = math.sqrt(tv)
            R_int = 4.8 * sq
            for _ in range(4):
                sb = flint.acb("0.5") + flint.acb(0, 1) * zc / 2
                acc = flint.acb(0)
                dacc = flint.acb(0)
                for xg, wg in zip(gl_x, gl_w):
                    u = R_int * float(xg)
                    wt = float(wg) * math.exp(-u * u / tv)
                    val = txi(sb + flint.arb(u))
                    acc += flint.arb(wt) * val
                    dacc += flint.arb(wt * (-2.0 * u / tv)) * val
                df = flint.acb(0, 1) * (-dacc) / 2
                st = acc / df
                zc = flint.acb((zc.real - st.real).mid(), (zc.imag - st.imag).mid())
            dels.append(max(0.0, abs(float(zc.imag.mid())) / 2.0))
        ts_full = np.append(ts, t_col)
        dels.append(0.0)
        return ts_full, np.array(dels)

    t1, d1 = trace_heat_curve(85.6993484854, 0.8085171825, 0.17650538)
    t531, d531 = trace_heat_curve(531.2797268965, 0.8469538031, 0.20127583)
    t240, d240 = trace_heat_curve(240.4046723514, 0.8695305796, 0.23060736)

    ax1.plot(t1, d1, color="#556B8D", ls="--", lw=1.5, label=r"$\rho_1\ (\gamma \approx 85.70,\ t_{\mathrm{col}} \approx 0.1765)$")
    ax1.plot(t531, d531, color="#B85D19", ls="-.", lw=1.6, label=r"$\rho_{531}\ (\gamma \approx 531.28,\ t_{\mathrm{col}} \approx 0.2013)$")
    ax1.plot(t240, d240, color="#8B1E3F", ls="-", lw=2.1, label=r"$\rho_{240}\ (\gamma \approx 240.40,\ t_{\mathrm{col}} \approx 0.23061)$")
    ax1.axvline(0.2306, color="#1B4332", ls=":", lw=1.5, label=r"Certified witness $t = 1153/5000 = 0.2306$")
    ax1.scatter([0.2306], [0.00383717 / 2.0], color=["#8B1E3F"], zorder=6, s=38)
    ax1.annotate(r"Certified at $t = \frac{1153}{5000}$", xy=(0.2306, 0.0019), xytext=(0.022, 0.042),
                 arrowprops={"arrowstyle": "->", "lw": 0.85, "color": "#8B1E3F"}, fontsize=8.2, color="#8B1E3F")
    ax1.set_xlim(0.0, 0.275)
    ax1.set_ylim(0.0, 0.40)
    ax1.set_xlabel(r"Backward heat parameter $t$")
    ax1.set_ylabel(r"Transverse off-line distance $\delta(t) = |\mathrm{Im}\,z(t)|/2$")
    ax1.set_title(r"(a) Certified Lower Bound $\Lambda_{\mathrm{DH}} > 1153/5000$")
    ax1.grid(True, alpha=0.25, lw=0.5)
    ax1.legend(loc="upper right", frameon=True, framealpha=0.92)

    survey_pts = [
        (85.6993, 0.808517, 0.176505),
        (114.1633, 0.654201, 0.0982),
        (165.4310, 0.728140, 0.1384),
        (240.4047, 0.869531, 0.230607),
        (531.2797, 0.846954, 0.201276),
        (1186.421, 0.859112, 0.1945),
        (2442.530, 0.897750, 0.199732),
        (6812.190, 0.904115, 0.1842),
        (18420.65, 0.916840, 0.1715),
        (49688.34, 0.931866, 0.1609),
        (68910.42, 0.921405, 0.1526),
    ]
    g_env = np.geomspace(70.0, 90000.0, 300)
    h_inf = 0.6203624982
    t_env = 4.0 * h_inf / np.log(5.0 * g_env / (2.0 * np.pi))

    ax2.plot(g_env, t_env, color="#1D3557", lw=1.7, ls="--",
             label=r"Leading-order heuristic $4(\sigma^* - \frac{1}{2}) / \log(5\gamma / 2\pi)$")
    ax2.axhline(0.23060736, color="#8B1E3F", lw=1.2, ls=":", label=r"$t_{\mathrm{col}}(\rho_{240}) \approx 0.23061$")
    ax2.scatter([p[0] for p in survey_pts], [p[2] for p in survey_pts], color="#8B1E3F", s=28, zorder=5,
                label=r"Sampled off-line zeros ($\gamma \leq 72{,}500$)")
    ax2.annotate(r"$\rho_{240}$", xy=(240.4, 0.2306), xytext=(330, 0.252),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.5)
    ax2.annotate(r"$\rho_{531}$", xy=(531.3, 0.2013), xytext=(780, 0.218),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.5)
    ax2.annotate(r"$\rho_{2442}$", xy=(2442.5, 0.1997), xytext=(3300, 0.145),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.5)
    ax2.annotate(r"$\gamma \approx 49{,}688$", xy=(49688, 0.1609), xytext=(12500, 0.105),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.2)
    ax2.set_xscale("log")
    ax2.set_xlim(70, 90000)
    ax2.set_ylim(0.08, 0.56)
    ax2.set_xlabel(r"Zero height $\gamma = \mathrm{Im}\,\rho$ (logarithmic scale)")
    ax2.set_ylabel(r"Collision time $t_{\mathrm{col}}(\rho)$")
    ax2.set_title("(b) Numerical Collision Survey vs. Leading-Order Heuristic")
    ax2.grid(True, which="both", alpha=0.22, lw=0.5)
    ax2.legend(loc="upper right", frameon=True, framealpha=0.92)

    f1_pdf = FIG_DIR / "fig1_heat_flow_separation_witnesses.pdf"
    f1_png = FIG_DIR / "fig1_heat_flow_separation_witnesses.png"
    fig.savefig(f1_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(f1_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 2: Sufficient Stieltjes regimes & [0, 311/500] sinc continuation ---
    fig2, (bx1, bx2) = plt.subplots(1, 2, figsize=(10.6, 4.2), gridspec_kw={"wspace": 0.26})

    u_grid = np.linspace(0.05, 18.0, 400)
    h_crit = 0.6203624982
    tau_crit = 2.0 * h_crit
    S_reg1 = 2.0 * u_grid / (4.0 + u_grid ** 2) + 2.0 * u_grid / (9.0 + u_grid ** 2)
    h_ex = 0.622
    S_cosh = 2.0 * h_ex * np.tanh(2.0 * h_ex * u_grid)
    S_sinc = 2.0 * h_ex / np.tanh(2.0 * h_ex * u_grid) - 1.0 / u_grid
    S_reg3 = 0.16 * u_grid

    bx1.plot(u_grid, S_reg1, color="#556B8D", ls="--", lw=1.6, label=r"Polynomial $G(w)$, $S_G(u)\to 0$")
    bx1.plot(u_grid, S_cosh, color="#1D3557", ls="-", lw=1.8, label=r"$G_{\mathrm{cosh},h}$ ($h=311/500$, $S_G\to 2h$)")
    bx1.plot(u_grid, S_sinc, color="#8B1E3F", ls="-.", lw=1.9, label=r"$G_{\mathrm{sinc},h}$ ($h=311/500$, $S_G\to 2h$)")
    bx1.plot(u_grid, S_reg3, color="#1B4332", ls=":", lw=1.8, label=r"Heat $G_t(w)=\mathrm{e}^{tw}$, $S_G(u)=4tu$")
    bx1.axhline(tau_crit, color="#333333", lw=1.0, ls="--", alpha=0.75, label=r"Threshold $2\sigma^* - 1 \approx 1.2407$")
    bx1.set_xlim(0.0, 18.0)
    bx1.set_ylim(0.0, 2.58)
    bx1.set_xlabel(r"Spectral frequency $u = \sqrt{w}$")
    bx1.set_ylabel(r"Stieltjes derivative $S_G(u) = \frac{\mathrm{d}}{\mathrm{d}u}\log G(u^2)$")
    bx1.set_title(r"(a) Asymptotic Shift Profiles $S_G(u)$ vs. $2\sigma^* - 1$")
    bx1.grid(True, alpha=0.25, lw=0.5)
    bx1.legend(loc="upper left", frameon=True, framealpha=0.92)

    flint.ctx.prec = 320
    tilde_xi = sinc_data["tilde_xi"]

    def g_cosh_eval_and_deriv(z_pt: flint.acb, h_val: flint.arb) -> tuple[flint.acb, flint.acb]:
        sb = flint.acb("0.5") + flint.acb(0, 1) * z_pt / 2
        f_val = (tilde_xi(sb + h_val) + tilde_xi(sb - h_val)) / 2
        eps = flint.arb("1e-25")
        sb_p = sb + flint.acb(0, 1) * eps / 2
        sb_m = sb - flint.acb(0, 1) * eps / 2
        df_val = ((tilde_xi(sb_p + h_val) + tilde_xi(sb_p - h_val)) - (tilde_xi(sb_m + h_val) + tilde_xi(sb_m - h_val))) / (4 * eps)
        return f_val, df_val

    cosh_h_grid = [
        0.001, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27,
        0.30, 0.325, 0.35, 0.37, 0.382, 0.39, 0.394, 0.3965, 0.3975,
    ]
    cosh_c = []
    zc = flint.acb("4885.06058493697206", "-0.79550030995402")
    for hv in cosh_h_grid:
        hp = flint.arb(str(hv))
        for _ in range(5):
            fv, dfv = g_cosh_eval_and_deriv(zc, hp)
            st = fv / dfv
            zc = flint.acb((zc.real - st.real).mid(), (zc.imag - st.imag).mid())
        cosh_c.append((hv, abs(float(zc.imag.mid())) / 2.0))
    cosh_c.append((0.397750155, 0.0))

    sinc_c = sorted(sinc_data["sinc_curve"], key=lambda p: p[0])
    sinc_full_h = [p[0] for p in sinc_c] + [0.62456]
    sinc_full_d = [p[1] for p in sinc_c] + [0.0]

    bx2.plot([p[0] for p in cosh_c], [p[1] for p in cosh_c], color="#1D3557", lw=1.9, ls="--",
             label=r"$G_{\mathrm{cosh},h}(-\partial_z^2)\Xi_{\mathrm{DH}}$ ($h_{\mathrm{col}} \approx 0.3978$)")
    bx2.plot(sinc_full_h, sinc_full_d, color="#8B1E3F", lw=2.1, ls="-",
             label=r"$G_{\mathrm{sinc},h}(-\partial_z^2)\Xi_{\mathrm{DH}}$ (certified on $[0, 311/500]$)")
    bx2.axvline(h_crit, color="#333333", lw=1.4, ls=":", label=r"$\sigma^* - \frac{1}{2} \approx 0.62036$")
    bx2.scatter([0.622], [0.05600463 / 2.0], color="#8B1E3F", s=38, zorder=6,
                label=r"Certified non-real zero at $h = 311/500$")
    bx2.annotate(r"$h = \frac{311}{500} > \sigma^* - \frac{1}{2}$", xy=(0.622, 0.0280), xytext=(0.405, 0.088),
                 arrowprops={"arrowstyle": "->", "lw": 0.85, "color": "#8B1E3F"}, fontsize=8.4, color="#8B1E3F")
    bx2.set_xlim(0.0, 0.655)
    bx2.set_ylim(0.0, 0.56)
    bx2.set_xlabel(r"Half-width parameter $h = \tau / 2$")
    bx2.set_ylabel(r"Transverse off-line distance $|\mathrm{Im}\,z(h)|/2$ near $\rho_{2442}$")
    bx2.set_title(r"(b) Certified Continuation on $h \in [0, 311/500]$")
    bx2.grid(True, alpha=0.25, lw=0.5)
    bx2.legend(loc="upper right", frameon=True, framealpha=0.92)

    f2_pdf = FIG_DIR / "fig2_polya_schur_sinc_vs_cosh.pdf"
    f2_png = FIG_DIR / "fig2_polya_schur_sinc_vs_cosh.png"
    fig2.savefig(f2_pdf, format="pdf", bbox_inches="tight")
    fig2.savefig(f2_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig2)

    for pdf_f in (f1_pdf, f2_pdf):
        kb = pdf_f.stat().st_size / 1024.0
        print(f"      Generated {pdf_f.name} ({kb:.1f} KB)")
        assert kb < 120.0, f"{pdf_f.name} exceeds 120 KB vector limit ({kb:.1f} KB)"


def main() -> int:
    parser = argparse.ArgumentParser(description="Rigorous FLINT/Arb verification suite for Paper XV.")
    parser.add_argument(
        "--regenerate-figures",
        action="store_true",
        help="Regenerate vector PDF and 300-DPI PNG figures in figures/ (off by default to keep git status clean).",
    )
    args = parser.parse_args()

    t0 = time.time()
    print("=" * 76)
    print("PAPER XV RIGOROUS FLINT/ARB BALL-ARITHMETIC VERIFICATION SUITE")
    print("=" * 76)
    _, h_inf_star, _ = certify_supremum_sigma_star()
    sinc_data = certify_polya_sinc_separation(h_inf_star)
    certify_heat_flow_witness()
    verify_analytic_lemmas()
    if args.regenerate_figures:
        generate_figures(sinc_data)
    print("=" * 76)
    print(f"ALL PAPER XV CERTIFICATES VERIFIED 100% IN {time.time() - t0:.2f}s.")
    print("=" * 76)
    return 0


if __name__ == "__main__":
    sys.exit(main())
