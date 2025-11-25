import numpy as np
from scipy.fft import rfft, irfft
# from scipy.linalg import solve_banded
# from scipy.integrate import cumulative_trapezoid
# from dataclasses import dataclass
from .zdiff import diff1_z_k, diff2_z_k, diff3_z_k, dealias_mask_rfft


def nonlinear_terms(wk, dz, k_vec, dealias=True):
    wk = np.asarray(wk)
    Nz, Nk = wk.shape
    Nx = 2*(Nk-1)

    # z-derivs (in k-space)
    dwk_dz   = diff1_z_k(wk, dz)
    d2wk_dz2 = diff2_z_k(wk, dz)
    d3wk_dz3 = diff3_z_k(wk, dz)

    # k helpers (row-shaped)
    ik      = (1j * k_vec)[None, :]
    k2_row  = (k_vec**2)[None, :]
    inv_ik  = np.zeros_like(ik, dtype=complex)
    nz = k_vec != 0
    inv_ik[:, nz] = 1.0 / (1j * k_vec[nz])

    # compose in x-physical space
    P  = irfft(inv_ik * dwk_dz,    n=Nx, axis=1, workers=-1).real
    D1 = irfft(d2wk_dz2,           n=Nx, axis=1, workers=-1).real
    Q  = irfft(wk,                 n=Nx, axis=1, workers=-1).real
    D2 = irfft(inv_ik * d3wk_dz3,  n=Nx, axis=1, workers=-1).real
    D3 = irfft(-k2_row * wk,       n=Nx, axis=1, workers=-1).real
    D4 = irfft(ik * dwk_dz,        n=Nx, axis=1, workers=-1).real

    NL_phys = P*(D1 + D3) - Q*(D2 + D4)
    NLk = rfft(NL_phys, axis=1, workers=-1)

    if dealias:
        mask = dealias_mask_rfft(Nk)
        NLk[:, ~mask] = 0.0

    return NLk




