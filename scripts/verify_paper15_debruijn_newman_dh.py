#!/usr/bin/env python3
"""
Rigorous C-Arb Ball-Arithmetic Verification Suite for Paper XV (v15.2.0)
========================================================================
Manuscript: Asymptotic Shift Thresholds for Polya-Schur Multipliers and a
            Certified Lower Bound on the de Bruijn-Newman Constant of the
            Davenport-Heilbronn Function
Author:     Thijs Missiaen (ORCID: 0009-0007-3678-0810)

Default usage (verifies all FLINT/Arb certificates without touching figures/):
    python scripts/verify_paper15_debruijn_newman_dh.py

Optional figure regeneration (writes PDF/PNG plots to figures/):
    python scripts/verify_paper15_debruijn_newman_dh.py --regenerate-figures

Certifies in 320-bit, 384-bit, and 512-bit FLINT/Arb ball arithmetic (python-flint):
  [1/4] IVT enclosure of sigma* = sup{Re(rho) : D_DH(rho) = 0} (radius 10^-30),
        h_infty* = sigma* - 1/2, 2*sigma* - 1, 2*(sigma* - 1/2)^2, and Euler-product
        phase sign -(1 - i*kappa)/(1 + i*kappa) = exp(+2i*(pi/2 - arctan(kappa))).
  [2/4] Singularity-free Taylor-series + 256-arc Cauchy-circle Interval-Newton
        continuation across the ENTIRE compact interval h in [0, 5/8] (360
        overlapping h-balls including h = 0, plus radius 10^-30 enclosures at
        h = 5/8 = 0.625 and h = 625121/10^6 = 0.625121 near rho_44838, and
        critical-boundary enclosures at h_* = sigma* - 1/2 near both rho_2442
        and rho_44838), proving 5/8 < h_global*(G_sinc; Xi_DH) <= 2*sigma* - 1
        independently of Section 2 and N_off(G_{sinc,sigma*-1/2}(-d_z^2)Xi_DH) >= 8.
  [3/4] 512-bit Interval-Newton off-real zero enclosure for H_t^DH(z) at exact
        rational t = 28825920237 / 125000000000 = 0.230607361896 near rho_240
        (radius 10^-32), proving Lambda_DH > 28825920237 / 125000000000.
  [4/4] Explicit analytical constants for Lemma 4.2, Proposition 4.3, Proposition 3.3
        (branch-locking & Mellin residue constants on Re(w) >= 1.07), Remark 4.5
        (Dobner Gaussian l^1-dominance threshold), and Proposition 2.7 (tower-lacunary).
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
    sinc_upper = 2 * h_inf_star
    debruijn_bound = 2 * (h_inf_star ** 2)

    print(f"      sigma*           = {sigma_star.str(36, radius=True)}")
    print(f"      sigma* - 1/2     = {h_inf_star.str(36, radius=True)}")
    print(f"      2*sigma* - 1     = {sinc_upper.str(36, radius=True)}")
    print(f"      2*(sigma*-1/2)^2 = {debruijn_bound.str(36, radius=True)}")
    return sigma_star, h_inf_star, debruijn_bound


# ============================================================================
# 2. SINGULARITY-FREE CONTINUUM CERTIFICATE ON ALL h in [0, 5/8] (rho_44838)
# ============================================================================
def certify_polya_sinc_separation(h_inf_star: flint.arb) -> dict[str, object]:
    print("\n[2/4] Certifying singularity-free Interval-Newton cover across ALL h in [0, 5/8] (rho_44838)...")
    flint.ctx.prec = 384
    sqrt5 = flint.arb(5).sqrt()
    kap = ((flint.arb(10) - 2 * sqrt5).sqrt() - 2) / (sqrt5 - 1)
    pi_a = flint.arb.pi()
    h_58 = flint.arb(5) / flint.arb(8)
    h_sharp = flint.arb(625121) / flint.arb(1000000)
    assert bool(h_58 > h_inf_star and h_sharp > h_58)

    s_ref = flint.acb(flint.arb("0.5"), flint.arb("44838.46"))

    def tilde_xi_44838(s: flint.acb) -> flint.acb:
        lp = (flint.acb(5) / pi_a).log()
        lg = ((s + 1) / 2).lgamma() - ((s_ref + 1) / 2).lgamma()
        rg = (((s - s_ref) / 2) * lp + lg).exp()
        fv = ((-s * flint.acb(5).log()).exp()) * (
            s.zeta(flint.acb(1) / 5) + kap * s.zeta(flint.acb(2) / 5)
            - kap * s.zeta(flint.acb(3) / 5) - s.zeta(flint.acb(4) / 5)
        )
        return rg * fv

    # Certify M_3 := sup_{|w - s_ref| = 3} |tilde_xi(w)| via 256 overlapping interval arcs
    M_3 = flint.arb(0)
    half_step = flint.arb("0.5001")
    for k in range(256):
        th = 2 * pi_a * flint.arb(flint.arb(k) + flint.arb("0.5"), half_step) / 256
        w = flint.arb(3) * (flint.acb(0, 1) * th).exp()
        M_3 = M_3.max(abs(tilde_xi_44838(s_ref + w)).upper())
    print(f"      256-arc Cauchy bound on |w - s_ref| = 3: M_3 <= {float(M_3.upper()):.4e}")

    ORD_PT = 128
    ORD_SL = 52
    flint.ctx.cap = ORD_PT
    ss = flint.acb_series([s_ref, 1])
    lp = (flint.acb(5) / pi_a).log()
    lg = ((ss + 1) / 2).lgamma() - ((s_ref + 1) / 2).lgamma()
    rg = (((ss - s_ref) / 2) * lp + lg).exp()
    fv = ((-ss * flint.acb(5).log()).exp()) * (
        ss.zeta(flint.acb(1) / 5) + kap * ss.zeta(flint.acb(2) / 5)
        - kap * ss.zeta(flint.acb(3) / 5) - ss.zeta(flint.acb(4) / 5)
    )
    jet = rg * fv
    cks = [jet[k] for k in range(ORD_PT)]
    int_cks = [cks[k] / (k + 1) for k in range(ORD_PT)]
    binom_even = [[flint.arb(math.comb(k, 2 * m)) / (2 * m + 1) for m in range(k // 2 + 1)] for k in range(ORD_SL)]
    binom_dh = [[flint.arb(math.comb(k, 2 * m) * 2 * m) / (2 * m + 1) for m in range(1, k // 2 + 1)] for k in range(ORD_SL)]

    def eval_sinc_slice(ds: flint.acb, h: flint.arb) -> tuple[flint.acb, flint.acb, flint.acb]:
        # Singularity-free power series on all h in R (including h = 0, zero division by h)
        rm = abs(ds).upper() + abs(h).upper()
        assert bool(rm < flint.arb("1.2")), "Trajectory left inner disk |w - s_ref| < 1.2"
        q = rm / 3
        err0 = M_3 * (q ** ORD_SL) / (1 - q)
        err1 = (M_3 / 3) * ORD_SL * (q ** (ORD_SL - 1)) / ((1 - q) ** 2)
        h2 = h * h
        h2p = [flint.arb(1)]
        for _ in range(ORD_SL // 2):
            h2p.append(h2p[-1] * h2)
        dsp = [flint.acb(1)]
        for _ in range(ORD_SL):
            dsp.append(dsp[-1] * ds)
        v0 = flint.acb(0)
        v1 = flint.acb(0)
        vh = flint.acb(0)
        for k in range(ORD_SL):
            ck = cks[k]
            s0 = flint.acb(0)
            for m, bm in enumerate(binom_even[k]):
                s0 += bm * dsp[k - 2 * m] * h2p[m]
            v0 += ck * s0
            if k >= 1:
                s1 = flint.acb(0)
                for m, bm in enumerate(binom_even[k - 1]):
                    s1 += bm * dsp[k - 1 - 2 * m] * h2p[m]
                v1 += (ck * k) * s1
            if k >= 2:
                sh = flint.acb(0)
                for idx, bm in enumerate(binom_dh[k]):
                    m = idx + 1
                    sh += bm * dsp[k - 2 * m] * (h * h2p[m - 1])
                vh += ck * sh
        return (
            v0 + flint.acb(flint.arb(0, err0.upper()), flint.arb(0, err0.upper())),
            (flint.acb(0, 1) / 2) * (v1 + flint.acb(flint.arb(0, err1.upper()), flint.arb(0, err1.upper()))),
            vh + flint.acb(flint.arb(0, err1.upper()), flint.arb(0, err1.upper())),
        )

    def eval_xi_ball(u: flint.acb) -> flint.acb:
        acc = flint.acb(0)
        for k in range(ORD_PT - 1, -1, -1):
            acc = acc * u + cks[k]
        ru = abs(u).upper()
        assert bool(ru < flint.arb("1.2"))
        q = ru / 3
        err = M_3 * (q ** ORD_PT) / (1 - q)
        eb = flint.arb(0, err.upper())
        return acc + flint.acb(eb, eb)

    def eval_int_diff(ds: flint.acb, h: flint.arb) -> flint.acb:
        up = ds + h
        um = ds - h
        ap = flint.acb(0)
        am = flint.acb(0)
        for k in range(ORD_PT - 1, -1, -1):
            ap = ap * up + int_cks[k]
            am = am * um + int_cks[k]
        r_seg = abs(up).upper().max(abs(um).upper())
        assert bool(r_seg < flint.arb("1.2"))
        q = r_seg / 3
        err = M_3 * (q ** ORD_PT) / (1 - q)
        eb = flint.arb(0, err.upper())
        return (ap * up - am * um) / (2 * h) + flint.acb(eb, eb)

    segments = [
        (flint.arb(0), flint.arb("0.55"), 120, flint.arb("8e-3")),
        (flint.arb("0.55"), flint.arb("0.620"), 100, flint.arb("3e-3")),
        (flint.arb("0.620"), flint.arb("0.6247"), 80, flint.arb("8e-4")),
        (flint.arb("0.6247"), h_58, 60, flint.arb("2e-4")),
    ]
    z_tr = flint.acb("89676.9504213421348454833650904697546", "-0.83227516635938639129247749946982753")
    sinc_curve_44838: list[tuple[float, float]] = []
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
                sinc_curve_44838.append((0.0, abs(float(z_tr.imag.mid())) / 2.0))
            if seg_idx == len(segments) - 1 and k_s == n_balls - 1:
                assert H_k.contains(h_58), "Final h-ball must contain h = 5/8"

            for _ in range(3):
                ds0 = flint.acb("0.5") + flint.acb(0, 1) * z_tr / 2 - s_ref
                f0 = eval_int_diff(ds0, h_mid)
                df0 = flint.acb(0, 1) * (eval_xi_ball(ds0 + h_mid) - eval_xi_ball(ds0 - h_mid)) / (4 * h_mid)
                st = f0 / df0
                z_tr = flint.acb((z_tr.real - st.real).mid(), (z_tr.imag - st.imag).mid())

            ds_mid = flint.acb("0.5") + flint.acb(0, 1) * z_tr / 2 - s_ref
            B_k = z_tr + flint.acb(flint.arb(0, r_box), flint.arb(0, r_box))
            ds_B = flint.acb("0.5") + flint.acb(0, 1) * B_k / 2 - s_ref
            g_mid, _, _ = eval_sinc_slice(ds_mid, h_mid)
            _, _, dh_g = eval_sinc_slice(ds_mid, H_k)
            _, dz_g_Bk, _ = eval_sinc_slice(ds_B, H_k)

            N_Bk = z_tr - (g_mid + (H_k - h_mid) * dh_g) / dz_g_Bk
            assert B_k.contains_interior(N_Bk) and bool(N_Bk.imag < flint.arb(0)), (
                f"Continuum sinc cert failed at segment {seg_idx}, slice {k_s}"
            )
            if seg_idx == 0 and k_s == 0:
                z_rho0 = flint.acb("89676.9504213421348454833650904697546", "-0.83227516635938639129247749946982753")
                B_rho0 = z_rho0 + flint.acb(flint.arb(0, flint.arb("1e-12")), flint.arb(0, flint.arb("1e-12")))
                ds_r0 = flint.acb("0.5") + flint.acb(0, 1) * z_rho0 / 2 - s_ref
                ds_Br0 = flint.acb("0.5") + flint.acb(0, 1) * B_rho0 / 2 - s_ref
                _, dz_r0, _ = eval_sinc_slice(ds_Br0, flint.arb(0))
                N_rho0 = z_rho0 - eval_xi_ball(ds_r0) / dz_r0
                assert B_rho0.contains_interior(N_rho0) and B_k.contains(B_rho0), (
                    "Initial h = 0 zero enclosure B_rho0 (rho_44838) failed or not inside B_1"
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
            if k_s % 8 == 0 or seg_idx >= 2:
                sinc_curve_44838.append((float(h_mid.mid()), abs(float(z_tr.imag.mid())) / 2.0))

    def certify_pt_44838(zc_pt: flint.acb, h_val: flint.arb, r_pt: flint.arb) -> flint.acb:
        Bp = zc_pt + flint.acb(flint.arb(0, r_pt), flint.arb(0, r_pt))
        dsp = flint.acb("0.5") + flint.acb(0, 1) * zc_pt / 2 - s_ref
        dsB = flint.acb("0.5") + flint.acb(0, 1) * Bp / 2 - s_ref
        fp = eval_int_diff(dsp, h_val)
        dfp = flint.acb(0, 1) * (eval_xi_ball(dsB + h_val) - eval_xi_ball(dsB - h_val)) / (4 * h_val)
        Np = zc_pt - fp / dfp
        assert Bp.contains_interior(Np) and bool(Np.imag < flint.arb(0)), "Point cert failed for rho_44838"
        return Bp

    zc_58 = flint.acb(
        "89676.898678875343559173771577620034736",
        "-0.01056386075347512714311622077705919110",
    )
    B_58 = certify_pt_44838(zc_58, h_58, flint.arb("1e-30"))
    assert last_B_k is not None and last_B_k.contains(B_58), "B_58 must lie inside B_360 for branch continuity"
    sinc_curve_44838.append((0.625, abs(float(zc_58.imag.mid())) / 2.0))

    zc_sharp = flint.acb(
        "89676.898662653368485342445005525886056",
        "-0.00055069865930981622151740114454506183",
    )
    certify_pt_44838(zc_sharp, h_sharp, flint.arb("1e-30"))
    sinc_curve_44838.append((0.625121, abs(float(zc_sharp.imag.mid())) / 2.0))

    zc_44838_hstar = flint.acb(
        "89676.899296681820873729194536307475930",
        "-0.06634905818151495877633209560664307483",
    )
    certify_pt_44838(zc_44838_hstar, h_inf_star, flint.arb("1e-27"))

    # Also certify rho_2442 at h_* = sigma* - 1/2 and h = 311/500 (320-bit, 128-arc Cauchy bound on |w - s_ref2| = 2):
    flint.ctx.prec = 320
    s_ref2 = flint.acb(flint.arb("0.71"), flint.arb("2442.49"))
    sqrt5_320 = flint.arb(5).sqrt()
    kap_320 = ((flint.arb(10) - 2 * sqrt5_320).sqrt() - 2) / (sqrt5_320 - 1)
    pi_320 = flint.arb.pi()

    def tilde_xi_2442(s: flint.acb) -> flint.acb:
        lp2 = (flint.acb(5) / pi_320).log()
        lg2 = ((s + 1) / 2).lgamma() - ((s_ref2 + 1) / 2).lgamma()
        rg2 = (((s - s_ref2) / 2) * lp2 + lg2).exp()
        fv2 = ((-s * flint.acb(5).log()).exp()) * (
            s.zeta(flint.acb(1) / 5) + kap_320 * s.zeta(flint.acb(2) / 5)
            - kap_320 * s.zeta(flint.acb(3) / 5) - s.zeta(flint.acb(4) / 5)
        )
        return rg2 * fv2

    M2_2442 = flint.arb(0)
    for k in range(128):
        th2 = 2 * pi_320 * flint.arb(flint.arb(k) + flint.arb("0.5"), flint.arb("0.5001")) / 128
        w2 = flint.arb(2) * (flint.acb(0, 1) * th2).exp()
        M2_2442 = M2_2442.max(abs(tilde_xi_2442(s_ref2 + w2)).upper())

    def certify_pt_2442(zc2: flint.acb, h2_val: flint.arb, r2_pt: flint.arb) -> None:
        flint.ctx.cap = 96
        sb2 = flint.acb("0.5") + flint.acb(0, 1) * zc2 / 2
        ss2 = flint.acb_series([sb2, 1])
        lp2 = (flint.acb(5) / pi_320).log()
        lg2 = ((ss2 + 1) / 2).lgamma() - ((s_ref2 + 1) / 2).lgamma()
        rg2 = (((ss2 - s_ref2) / 2) * lp2 + lg2).exp()
        fv2 = ((-ss2 * flint.acb(5).log()).exp()) * (
            ss2.zeta(flint.acb(1) / 5) + kap_320 * ss2.zeta(flint.acb(2) / 5)
            - kap_320 * ss2.zeta(flint.acb(3) / 5) - ss2.zeta(flint.acb(4) / 5)
        )
        j2 = rg2 * fv2
        R_av2 = flint.arb(2) - abs(sb2 - s_ref2).upper()
        q2 = abs(h2_val).upper() / R_av2
        err2 = M2_2442 * (q2 ** 96) / (1 - q2)
        f2 = sum((j2[2 * m] / (2 * m + 1)) * (h2_val ** (2 * m)) for m in range(48)) + flint.acb(
            flint.arb(0, err2.upper()), flint.arb(0, err2.upper())
        )
        B2 = zc2 + flint.acb(flint.arb(0, r2_pt), flint.arb(0, r2_pt))
        sb_B2 = flint.acb("0.5") + flint.acb(0, 1) * B2 / 2
        df2 = flint.acb(0, 1) * (tilde_xi_2442(sb_B2 + h2_val) - tilde_xi_2442(sb_B2 - h2_val)) / (4 * h2_val)
        N2 = zc2 - f2 / df2
        assert B2.contains_interior(N2) and bool(N2.imag < flint.arb(0)), "Point cert failed for rho_2442"

    zc_2442_hstar = flint.acb(
        "4884.88734363569522396877972427073875",
        "-0.071915260298399708621062650316073866",
    )
    certify_pt_2442(zc_2442_hstar, h_inf_star, flint.arb("1e-27"))

    zc_2442_622 = flint.acb(
        "4884.88665419175715643333300197615178",
        "-0.056004628860858072317639148896287841",
    )
    certify_pt_2442(zc_2442_622, flint.arb(311) / flint.arb(500), flint.arb("1e-30"))

    print(f"      Continuum & branch-overlap certified across ALL h in [0, 5/8] ({total_balls} balls, 359 joins): True")
    print(f"      z_44838(5/8)        = {zc_58.real.str(36, radius=False)} + ({zc_58.imag.str(36, radius=False)})j (r=10^-30)")
    print(f"      z_44838(625121/10^6)= {zc_sharp.real.str(36, radius=False)} + ({zc_sharp.imag.str(36, radius=False)})j (r=10^-30)")
    print(f"      Critical boundary h_* = sigma* - 1/2 certified for rho_2442 and rho_44838 (N_off >= 8, r=10^-27): True")

    return {
        "z_44838_58": zc_58,
        "sinc_curve_44838": sinc_curve_44838,
        "tilde_xi_2442": tilde_xi_2442,
    }


# ============================================================================
# 3. 512-BIT INTERVAL-NEWTON HEAT WITNESS AT t = 28825920237 / 125000000000
# ============================================================================
def certify_heat_flow_witness() -> flint.acb:
    print("\n[3/4] Certifying 512-bit Interval-Newton witness at t = 28825920237/125000000000 = 0.230607361896...")
    flint.ctx.prec = 512
    sqrt5 = flint.arb(5).sqrt()
    kap = ((flint.arb(10) - 2 * sqrt5).sqrt() - 2) / (sqrt5 - 1)
    pi_a = flint.arb.pi()
    pi5 = flint.acb(pi_a / 5)

    n_terms = 34
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

    U_top = flint.arb(23) / flint.arb(10)
    t_heat = flint.arb(28825920237) / flint.arb(125000000000)
    tol = flint.arb("1e-136")

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
        Im = 2 * flint.acb.integral(ker, flint.acb(0), flint.acb(U_top), abs_tol=tol, rel_tol=tol, eval_limit=300000)
        eb = tail_ball(abs(z.imag).upper(), order)
        rad = flint.arb(0, eb)
        return Im + flint.acb(rad, rad)

    zc = flint.acb(
        "480.6902310425074727072038413480701557071",
        "-9.535301365553738949062115339251329328331e-7",
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
    rem_d1 = (abs(dz).upper() ** 3) * M4 / 6
    df_B = dfz + dz * (d2fz + (dz / 2) * d3fz) + flint.acb(flint.arb(0, rem_d1.upper()), flint.arb(0, rem_d1.upper()))
    N_B = zc - fz / df_B

    ok_in = B.contains_interior(N_B)
    ok_im = bool(N_B.imag < flint.arb(0))
    assert ok_in and ok_im, "Interval-Newton failed for rho_240 at t = 0.230607361896"
    print(f"      z_240(0.230607361896) = {zc.real.str(36, radius=False)} + ({zc.imag.str(36, radius=False)})j")
    print(f"      N(B) c int(B) (r=10^-32): {ok_in},  Im(z) < 0: {ok_im}")
    return zc


# ============================================================================
# 4. EXPLICIT CONSTANTS FOR PROPOSITION 3.3, PROPOSITION 4.3 & PROPOSITION 2.7
# ============================================================================
def verify_analytic_lemmas(sigma_star: flint.arb) -> None:
    print("\n[4/4] Verifying explicit analytic bounds (Abel, Dirichlet tail, Prop. 3.3, Dobner l^1, lacunary)...")
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

    # 4c. Rigorous Arb constants for Proposition 3.3 (branch-locking & Mellin residue on Re(w) >= 1.07):
    flint.ctx.prec = 256
    s5_a = flint.arb(5).sqrt()
    kap_a = ((flint.arb(10) - 2 * s5_a).sqrt() - 2) / (s5_a - 1)
    pi_a = flint.arb.pi()
    theta_kap = pi_a - 2 * kap_a.atan()
    two_pi_minus_theta = 2 * pi_a - theta_kap

    sig_107 = flint.arb(107) / flint.arb(100)
    sc_107 = flint.acb(sig_107)
    zeta_107 = sc_107.zeta().real
    L_quad_107 = (
        (flint.acb(5) ** (-sc_107))
        * (sc_107.zeta(flint.acb(1) / 5) - sc_107.zeta(flint.acb(2) / 5) - sc_107.zeta(flint.acb(3) / 5) + sc_107.zeta(flint.acb(4) / 5))
    ).real
    # Using 2*arctan(x) < log((1+x)/(1-x)) for x in (0,1), sum_{p=2,3 mod 5} 2*arctan(p^-1.07) < log((1-5^-1.07)*zeta(1.07)/L(1.07,chi_5,quad)):
    phi_upper_107 = ((1 - flint.arb(5) ** (-sig_107)) * zeta_107 / L_quad_107).log()
    L_min_107 = (1 + flint.arb(5) ** (-sig_107)) * flint.acb(2 * sig_107).zeta().real / zeta_107
    c_star = (flint.arb(2) ** (1 - sigma_star) - 1) * (1 - sigma_star).gamma()
    assert bool(phi_upper_107 < flint.arb("3.33") and two_pi_minus_theta - phi_upper_107 > flint.arb("0.37"))
    assert bool(L_min_107 > flint.arb(3) / flint.arb(25) and c_star > flint.arb("0.71"))
    print(
        f"      Prop. 3.3 branch lock (sigma>=1.07): max|Phi|<={float(phi_upper_107.upper()):.4f} < 3.33 < "
        f"2*pi-theta_kap={float(two_pi_minus_theta.lower()):.4f}, |L|>={float(L_min_107.lower()):.4f}, c_*={float(c_star.mid()):.4f} (OK)"
    )

    # 4d. Gauss-Weierstrass half-line identity & Dobner Gaussian l^1-inflation threshold (Remark 4.5):
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
    # Dobner Riemann-Siegel-truncated l^1 tail sum_{2 <= n <= sqrt(5*gam/2*pi)} |a_DH(n)| * exp((t/4)*log^2 n) * n^(-1/2 - u_t):
    def dobner_rs_tail(t_v: float, gam_v: float) -> float:
        x_rs = 5.0 * gam_v / (2.0 * math.pi)
        n_rs = int(math.isqrt(int(x_rs)))
        ut = 0.25 * t_v * math.log(x_rs)
        tot = 0.0
        for n in range(2, n_rs + 1):
            r = n % 5
            if r == 0:
                continue
            an = 1.0 if r in (1, 4) else float(kap)
            ln = math.log(n)
            tot += an * math.exp(0.25 * t_v * ln * ln - (0.5 + ut) * ln)
        return tot
    assert dobner_rs_tail(0.50, 1.50e4) > 1.0 and dobner_rs_tail(0.50, 1.55e4) < 1.0
    assert dobner_rs_tail(0.25, 4.15e8) > 1.0 and dobner_rs_tail(0.25, 4.21e8) < 1.0
    print(f"      Gauss-Weierstrass half-line residual: {mp.nstr(fold_err, 8)} < 1e-35; Dobner RS l^1(0.25) ~ 4.18e8 (OK)")

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

    t_rat_plot = 0.230607361896
    ax1.plot(t1, d1, color="#556B8D", ls="--", lw=1.5, label=r"$\rho_1\ (\gamma \approx 85.70,\ t_{\mathrm{col}} \approx 0.1765)$")
    ax1.plot(t531, d531, color="#B85D19", ls="-.", lw=1.6, label=r"$\rho_{531}\ (\gamma \approx 531.28,\ t_{\mathrm{col}} \approx 0.2013)$")
    ax1.plot(t240, d240, color="#8B1E3F", ls="-", lw=2.1, label=r"$\rho_{240}\ (\gamma \approx 240.40,\ t_{\mathrm{col}} \approx 0.230607)$")
    ax1.axvline(t_rat_plot, color="#1B4332", ls=":", lw=1.5, label=r"Certified $t = 0.230607361896$")
    ax1.scatter([t_rat_plot], [9.5353e-7 / 2.0], color=["#8B1E3F"], zorder=6, s=38)
    ax1.annotate(r"Certified at $t = 0.230607361896$", xy=(t_rat_plot, 0.0005), xytext=(0.018, 0.042),
                 arrowprops={"arrowstyle": "->", "lw": 0.85, "color": "#8B1E3F"}, fontsize=8.0, color="#8B1E3F")
    ax1.set_xlim(0.0, 0.275)
    ax1.set_ylim(0.0, 0.40)
    ax1.set_xlabel(r"Backward heat parameter $t$")
    ax1.set_ylabel(r"Transverse off-line distance $\delta(t) = |\mathrm{Im}\,z(t)|/2$")
    ax1.set_title(r"(a) Certified Lower Bound $\Lambda_{\mathrm{DH}} > 0.230607361896$")
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
        (44838.48, 0.916138, 0.1592),
        (49688.34, 0.931866, 0.1609),
        (247535.1, 0.946120, 0.1418),
    ]
    g_env = np.geomspace(70.0, 350000.0, 320)
    h_inf = 0.6203624982
    t_env = 4.0 * h_inf / np.log(5.0 * g_env / (2.0 * np.pi))

    ax2.plot(g_env, t_env, color="#1D3557", lw=1.7, ls="--",
             label=r"Leading-order heuristic $4(\sigma^* - \frac{1}{2}) / \log(5\gamma / 2\pi)$")
    ax2.axhline(0.23060736, color="#8B1E3F", lw=1.2, ls=":", label=r"$t_{\mathrm{col}}(\rho_{240}) \approx 0.230607$")
    ax2.scatter([p[0] for p in survey_pts], [p[2] for p in survey_pts], color="#8B1E3F", s=26, zorder=5,
                label=r"Winding scan ($[0,10^4]$) & resonance peaks")
    ax2.annotate(r"$\rho_{240}$", xy=(240.4, 0.2306), xytext=(340, 0.256),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.4)
    ax2.annotate(r"$\rho_{531}$", xy=(531.3, 0.2013), xytext=(820, 0.218),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.4)
    ax2.annotate(r"$\rho_{2442}$", xy=(2442.5, 0.1997), xytext=(3100, 0.140),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.3)
    ax2.annotate(r"$\rho_{44838}$", xy=(44838.5, 0.1592), xytext=(11000, 0.104),
                 arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.2)
    ax2.set_xscale("log")
    ax2.set_xlim(70, 350000)
    ax2.set_ylim(0.08, 0.56)
    ax2.set_xlabel(r"Zero height $\gamma = \mathrm{Im}\,\rho$ (logarithmic scale)")
    ax2.set_ylabel(r"Collision time $t_{\mathrm{col}}(\rho)$")
    ax2.set_title("(b) Collision Survey vs. Leading-Order Heuristic")
    ax2.grid(True, which="both", alpha=0.22, lw=0.5)
    ax2.legend(loc="upper right", frameon=True, framealpha=0.92)

    f1_pdf = FIG_DIR / "fig1_heat_flow_separation_witnesses.pdf"
    f1_png = FIG_DIR / "fig1_heat_flow_separation_witnesses.png"
    fig.savefig(f1_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(f1_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 2: Sufficient Stieltjes regimes & [0, 5/8] sinc continuation ---
    fig2, (bx1, bx2) = plt.subplots(1, 2, figsize=(10.6, 4.2), gridspec_kw={"wspace": 0.26})

    u_grid = np.linspace(0.05, 18.0, 400)
    h_crit = 0.6203624982
    tau_crit = 2.0 * h_crit
    S_reg1 = 2.0 * u_grid / (4.0 + u_grid ** 2) + 2.0 * u_grid / (9.0 + u_grid ** 2)
    h_ex = 0.625
    S_cosh = 2.0 * h_ex * np.tanh(2.0 * h_ex * u_grid)
    S_sinc = 2.0 * h_ex / np.tanh(2.0 * h_ex * u_grid) - 1.0 / u_grid
    S_reg3 = 0.16 * u_grid

    bx1.plot(u_grid, S_reg1, color="#556B8D", ls="--", lw=1.6, label=r"Polynomial $G(w)$, $S_G(u)\to 0$")
    bx1.plot(u_grid, S_cosh, color="#1D3557", ls="-", lw=1.8, label=r"$G_{\mathrm{cosh},h}$ ($h=5/8$, $S_G\to 5/4$)")
    bx1.plot(u_grid, S_sinc, color="#8B1E3F", ls="-.", lw=1.9, label=r"$G_{\mathrm{sinc},h}$ ($h=5/8$, $S_G\to 5/4$)")
    bx1.plot(u_grid, S_reg3, color="#1B4332", ls=":", lw=1.8, label=r"Heat $G_t(w)=\mathrm{e}^{tw}$, $S_G(u)=4tu$")
    bx1.axhline(tau_crit, color="#333333", lw=1.0, ls="--", alpha=0.75, label=r"Threshold $2\sigma^* - 1 \approx 1.2407$")
    bx1.set_xlim(0.0, 18.0)
    bx1.set_ylim(0.0, 2.58)
    bx1.set_xlabel(r"Spectral frequency $u = \sqrt{w}$")
    bx1.set_ylabel(r"Stieltjes derivative $S_G(u) = \frac{\mathrm{d}}{\mathrm{d}u}\log G(u^2)$")
    bx1.set_title(r"(a) Asymptotic Shift Profiles $S_G(u)$ vs. $2\sigma^* - 1$")
    bx1.grid(True, alpha=0.25, lw=0.5)
    bx1.legend(loc="upper left", frameon=True, framealpha=0.92)

    flint.ctx.prec = 192
    tilde_xi_2442 = sinc_data["tilde_xi_2442"]

    def g_cosh_eval_and_deriv(z_pt: flint.acb, h_val: flint.arb) -> tuple[flint.acb, flint.acb]:
        sb = flint.acb("0.5") + flint.acb(0, 1) * z_pt / 2
        f_val = (tilde_xi_2442(sb + h_val) + tilde_xi_2442(sb - h_val)) / 2
        eps = flint.arb("1e-20")
        sb_p = sb + flint.acb(0, 1) * eps / 2
        sb_m = sb - flint.acb(0, 1) * eps / 2
        df_val = ((tilde_xi_2442(sb_p + h_val) + tilde_xi_2442(sb_p - h_val)) - (tilde_xi_2442(sb_m + h_val) + tilde_xi_2442(sb_m - h_val))) / (4 * eps)
        return f_val, df_val

    cosh_h_grid = [
        0.001, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27,
        0.30, 0.325, 0.35, 0.37, 0.382, 0.39, 0.394, 0.3965, 0.3975,
    ]
    cosh_c = [(0.0, 0.397750155)]
    zc = flint.acb("4885.06058493697206", "-0.79550030995402")
    for hv in cosh_h_grid:
        hp = flint.arb(str(hv))
        for _ in range(5):
            fv, dfv = g_cosh_eval_and_deriv(zc, hp)
            st = fv / dfv
            zc = flint.acb((zc.real - st.real).mid(), (zc.imag - st.imag).mid())
        cosh_c.append((hv, abs(float(zc.imag.mid())) / 2.0))
    cosh_c.append((0.397750155, 0.0))

    # Trace rho_2442 under G_sinc for comparison on Panel (b):
    sinc2442_grid = [0.05 * j for j in range(1, 12)] + [0.58, 0.60, 0.612, 0.618, 0.622, 0.624, 0.6244]
    sinc2442_c = [(0.0, 0.397750155)]
    zc2 = flint.acb("4885.06058493697206", "-0.79550030995402")
    gl24_x, gl24_w = np.polynomial.legendre.leggauss(24)
    for hv in sinc2442_grid:
        hp = flint.arb(str(hv))
        for _ in range(4):
            sb2 = flint.acb("0.5") + flint.acb(0, 1) * zc2 / 2
            fv = sum(flint.arb(float(wg) * 0.5) * tilde_xi_2442(sb2 + flint.arb(float(xg) * hv)) for xg, wg in zip(gl24_x, gl24_w))
            dfv = flint.acb(0, 1) * (tilde_xi_2442(sb2 + hp) - tilde_xi_2442(sb2 - hp)) / (4 * hp)
            st = fv / dfv
            zc2 = flint.acb((zc2.real - st.real).mid(), (zc2.imag - st.imag).mid())
        sinc2442_c.append((hv, abs(float(zc2.imag.mid())) / 2.0))
    sinc2442_c.append((0.624529, 0.0))

    sinc_c = sorted(sinc_data["sinc_curve_44838"], key=lambda p: p[0])
    sinc_full_h = [p[0] for p in sinc_c] + [0.6251213]
    sinc_full_d = [p[1] for p in sinc_c] + [0.0]

    bx2.plot([p[0] for p in cosh_c], [p[1] for p in cosh_c], color="#1D3557", lw=1.8, ls="--",
             label=r"$G_{\mathrm{cosh},h}$ on $\rho_{2442}$ ($h_{\mathrm{col}} \approx 0.3978$)")
    bx2.plot([p[0] for p in sinc2442_c], [p[1] for p in sinc2442_c], color="#B85D19", lw=1.6, ls="-.",
             label=r"$G_{\mathrm{sinc},h}$ on $\rho_{2442}$ ($h_{\mathrm{col}} \approx 0.62453$)")
    bx2.plot(sinc_full_h, sinc_full_d, color="#8B1E3F", lw=2.1, ls="-",
             label=r"$G_{\mathrm{sinc},h}$ on $\rho_{44838}$ (certified $[0, 5/8]$)")
    bx2.axvline(h_crit, color="#333333", lw=1.3, ls=":", label=r"$\sigma^* - \frac{1}{2} \approx 0.62036$")
    bx2.scatter([0.625], [0.01056386 / 2.0], color="#8B1E3F", s=36, zorder=6)
    bx2.annotate(r"Certified $h = \frac{5}{8}$", xy=(0.625, 0.0053), xytext=(0.415, 0.062),
                 arrowprops={"arrowstyle": "->", "lw": 0.85, "color": "#8B1E3F"}, fontsize=8.1, color="#8B1E3F")
    bx2.set_xlim(0.0, 0.662)
    bx2.set_ylim(0.0, 0.62)
    bx2.set_xlabel(r"Half-width parameter $h = \tau / 2$")
    bx2.set_ylabel(r"Transverse off-line distance $|\mathrm{Im}\,z(h)|/2$")
    bx2.set_title(r"(b) Certified Continuation Across $h \in [0, 5/8]$")
    bx2.grid(True, alpha=0.25, lw=0.5)
    bx2.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=7.8)

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
    print("PAPER XV RIGOROUS FLINT/ARB BALL-ARITHMETIC VERIFICATION SUITE (v15.2.0)")
    print("=" * 76)
    sigma_star, h_inf_star, _ = certify_supremum_sigma_star()
    sinc_data = certify_polya_sinc_separation(h_inf_star)
    certify_heat_flow_witness()
    verify_analytic_lemmas(sigma_star)
    if args.regenerate_figures:
        generate_figures(sinc_data)
    print("=" * 76)
    print(f"ALL PAPER XV CERTIFICATES VERIFIED 100% IN {time.time() - t0:.2f}s.")
    print("=" * 76)
    return 0


if __name__ == "__main__":
    sys.exit(main())
