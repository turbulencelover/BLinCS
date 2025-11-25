#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.fft import irfft, rfftfreq
from scipy.linalg import solve_banded
from scipy.integrate import cumulative_trapezoid

from .vel_base import BackgroundBase
from .forcing import ForcingBase, NonlinearForcing
from .nonlinear import nonlinear_terms


# ----------------- helper for one-sided BCs -----------------
def _apply_one_sided_BC(A, RHS, phi, dz, U0z):
    """
    Modify operator A and RHS in-place for one-sided boundary conditions:
      - bottom: w=0, dw/dz=0 (3rd order forward)
      - top: dw/dz = (phi/U0^2) w, d2w/dz2=0 (2nd order backward)
    """

    # lower BC
    A[:4, 0] = [0, 0, 1, -11]
    A[:3, 1] = [0, 0, 18]
    A[:2, 2] = [0, -9]
    A[0, 3]  = 2
    RHS[0] = 0
    RHS[1] = 0

    # upper BC
    A[-4:, -1] = [ 2, 3 - 2*phi/U0z[-1]**2*dz, 0, 0]
    A[-3:, -2] = [-5, -4, 0]
    A[-2:, -3] = [ 4,  1]
    A[-1, -4]  = -1
    RHS[-1] = 0
    RHS[-2] = 0


# ----------------- main solver -----------------
def simulate_BLwave(
    g_prime: float,
    N: float,
    background: BackgroundBase,
    forcing: ForcingBase,
    Nx: int = 1024,
    Nz: int = 1200,
    Lx: float = 100000,
    Lz: float = 500,
    Uwind: float = 20.0,
    nu_turb: float = 2.0,
    BC_method: str = "central",
    *,
    nonlinear_adv: bool = False,
    max_iter: int = 20,
    tol: float = 1e-6,
    dealias_NL: bool = True,
):
    """
    Linear solver with optional Picard iterations for nonlinear terms.
    """
    # --- grids ---
    dx = Lx / Nx
    x = np.linspace(-Lx/2, Lx/2, Nx, endpoint=False)
    z = np.linspace(1.0, Lz, Nz)
    dz = z[1] - z[0]
    xx, zz = np.meshgrid(x, z)

    # --- Fourier grid ---
    k_vec = 2 * np.pi * rfftfreq(Nx, dx)

    # --- background flow ---
    U0z, dU0z, d2U0z = background.eval(z)

    # --- forcing --- (universal call — works for linear AND nonlinear)
    Fxr, dfk_dz = forcing.generate(U0z, dU0z, k_vec, xx, zz, z)

    # --- allocate spectral solution ---
    Nk = len(k_vec)
    wk = np.zeros((Nz, Nk), dtype=complex)

    # --- base operators ---
    cD4 = np.zeros((5, Nz), dtype=complex)
    cD4[0, 2:] = 1
    cD4[1, 1:] = -4
    cD4[2, :]  = 6
    cD4[3, :-1] = -4
    cD4[4, :-2] = 1
    cD4[2, 0] = 7
    cD4 /= dz**4

    U0D2 = np.zeros((5, Nz), dtype=complex)
    U0D2[1, :] = U0z
    U0D2[2, :] = -2 * U0z
    U0D2[3, :] = U0z
    U0D2 /= dz**2

    nu_t = np.full((5, Nz), float(nu_turb))
    RHS = np.zeros(Nz, dtype=complex)
    
    # nonlinear forcing toggle (your preferred style)
    do_nl_forcing = (type(forcing) is NonlinearForcing)
    
    # nonlinear mode → force one-sided BCs
    if (nonlinear_adv or do_nl_forcing) and BC_method.lower() not in ("one-sided", "one-sided derivative", "sided derivative"):
        print("⚠️  Nonlinear mode requires one-sided BCs. Switching to one-sided.")
        BC_method = "one-sided"

    
   

    # ========== initial linear solve ==========
    for i, k in enumerate(k_vec):
        if i == 0:
            continue
        phi = g_prime + 1j * Uwind * N * np.sign(k)
        nuCoeff = 1j * nu_t / k

        RHS[:] = -dfk_dz[:, i]
        U0k2_d2U0 = np.zeros((5, Nz), dtype=complex)
        U0k2_d2U0[2, :] = -(U0z * k**2 + d2U0z)

        if BC_method.lower() in ("ghost-point", "ghost point", "central"):
            # Robin BC
            cD4[1, -1] = (-3 + phi*dz / (U0z[-1]**2 - phi*dz)) / dz**4
            cD4[2, -1] = ( 3 - 2*phi*dz / (U0z[-1]**2 - phi*dz)) / dz**4
            U0D2[2, -1] = (U0z[-1]*(-1 + phi*dz / (U0z[-1]**2 - phi*dz))) / dz**2
            A = U0D2 + U0k2_d2U0 + nuCoeff * cD4

        else:  # one-sided
            A = U0D2 + U0k2_d2U0 + nuCoeff * cD4
            _apply_one_sided_BC(A, RHS, phi, dz, U0z)

        wk[:, i] = solve_banded((2, 2), A, RHS)

    wk[:, 0] = 0
    

    # ========== nonlinear Picard iterations ==========
    if nonlinear_adv or do_nl_forcing:
        for it in range(max_iter):
            # recompute forcing if it depends on the state
            if do_nl_forcing:
                # reconstruct ur from current wk
                wr_it = np.real(irfft(wk, n=Nx, axis=1))
                dwr = np.zeros_like(wr_it)
                dwr[1:-1, :] = (wr_it[2:, :] - wr_it[:-2, :]) / (2 * dz)
                dwr[0, :]    =  wr_it[1, :] / (2 * dz)
                dwr[-1, :]   = (wr_it[-1, :] - wr_it[-2, :]) / (2 * dz)
                ur_it = -cumulative_trapezoid(dwr, x, axis=1, initial=0.0)
                
                Fxr, dfk_dz = forcing.generate(
                    U0z, dU0z, k_vec, xx, zz, z, state={"ur": ur_it} )
                
            # advective NL terms (can be off)
            NLk = nonlinear_terms(wk, dz, k_vec, dealias=dealias_NL) if nonlinear_adv else 0.0

            wk_new = np.zeros_like(wk)
            for i, k in enumerate(k_vec):
                if i == 0:
                    continue
                phi = g_prime + 1j * Uwind * N * np.sign(k)
                nuCoeff = 1j * nu_t / k

                RHS[:] = -dfk_dz[:, i] + (NLk[:, i] if isinstance(NLk, np.ndarray) else 0.0)
                U0k2_d2U0 = np.zeros((5, Nz), dtype=complex)
                U0k2_d2U0[2, :] = -(U0z * k**2 + d2U0z)
                A = U0D2 + U0k2_d2U0 + nuCoeff * cD4
                _apply_one_sided_BC(A, RHS, phi, dz, U0z)

                wk_new[:, i] = solve_banded((2, 2), A, RHS)

            wk_new[:, 0] = 0

            rel = np.linalg.norm(wk_new - wk) / (np.linalg.norm(wk) + 1e-14)
            print(f"[Picard] iter {it:02d}: rel change = {rel:.3e}")
            wk = wk_new
            if rel < tol:
                break
        else:
            print("⚠️ Picard did not converge within max_iter.")

    # ========== back to physical space ==========
    wr = np.real(irfft(wk, n=Nx, axis=1))

    dwr = np.zeros_like(wr)
    dwr[1:-1, :] = (wr[2:, :] - wr[:-2, :]) / (2 * dz)
    dwr[0, :]    =  wr[1, :] / (2 * dz)
    dwr[-1, :]   = (wr[-1, :] - wr[-2, :]) / (2 * dz)
    ur = -cumulative_trapezoid(dwr, x, axis=1, initial=0.0)

    eta = cumulative_trapezoid(wr[-1, :] / (U0z[-1] + ur[-1, :]), x, initial=0.0)

    return x, z, ur, wr, eta, Fxr, U0z
