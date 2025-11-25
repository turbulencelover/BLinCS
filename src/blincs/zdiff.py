#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 10:43:51 2025

@author: hkafiaba
"""

# zdiff.py
import numpy as np
from scipy.fft import rfft, irfft, rfftfreq  # optional, for users of dealias/test utils

# --------- Dealias mask (2/3 rule) ---------
def dealias_mask_rfft(Nk: int) -> np.ndarray:
    """
    Build 2/3-rule mask for an RFFT spectrum of length Nk = Nx//2 + 1.
    Keeps indices 0..floor((Nk-1)/3). True=keep, False=zero.
    """
    cutoff = (Nk - 1) // 3
    mask = np.zeros(Nk, dtype=bool)
    mask[:cutoff + 1] = True
    return mask

# --------- z-derivatives (2nd-order accurate) ---------
def diff1_z_k(wk: np.ndarray, dz: float) -> np.ndarray:
    """
    First derivative ∂/∂z acting on wk (Nz, Nk) where each column is a k-mode (RFFT in x).
    Uses centered 2nd-order in interior, one-sided 2nd-order at z-boundaries.
    """
    wk = np.asarray(wk)
    Nz, Nk = wk.shape
    d = np.empty_like(wk, dtype=wk.dtype)
    d[1:-1, :] = (wk[2:, :] - wk[:-2, :]) / (2.0 * dz)
    d[0,     :] = (-3.0*wk[0, :] + 4.0*wk[1, :] - wk[2, :]) / (2.0 * dz)
    d[-1,    :] = ( 3.0*wk[-1, :] - 4.0*wk[-2, :] + wk[-3, :]) / (2.0 * dz)
    return d

def diff2_z_k(wk: np.ndarray, dz: float) -> np.ndarray:
    """
    Second derivative ∂²/∂z² on wk (Nz, Nk). Centered 2nd-order interior; one-sided 2nd-order at boundaries.
    """
    wk = np.asarray(wk)
    Nz, Nk = wk.shape
    d2 = np.empty_like(wk, dtype=wk.dtype)
    d2[1:-1, :] = (wk[2:, :] - 2.0*wk[1:-1, :] + wk[:-2, :]) / (dz*dz)
    d2[0,  :]   = ( 2.0*wk[0, :] - 5.0*wk[1, :] + 4.0*wk[2, :] - wk[3, :] ) / (dz*dz)
    d2[-1, :]   = ( 2.0*wk[-1, :] - 5.0*wk[-2, :] + 4.0*wk[-3, :] - wk[-4, :] ) / (dz*dz)
    return d2

def diff3_z_k(wk: np.ndarray, dz: float) -> np.ndarray:
    """
    Third derivative ∂³/∂z³ on wk (Nz, Nk). 5-pt 2nd-order accurate.
    Requires Nz >= 6.
    - interior (i=2..Nz-3): (f_{i-2} - 2 f_{i-1} + 2 f_{i+1} - f_{i+2}) / (2 h^3)
    - forward/backward 5-pt at the two boundaries.
    """
    wk = np.asarray(wk)
    Nz, Nk = wk.shape
    if Nz < 6:
        raise ValueError("diff3_z_k requires Nz >= 6 for 5-point stencils.")
    h3 = 2.0 * dz**3
    d3 = np.empty_like(wk, dtype=wk.dtype)

    # interior
    d3[2:-2, :] = (wk[0:-4, :] - 2.0*wk[1:-3, :] + 2.0*wk[3:-1, :] - wk[4:, :]) / h3
    # forward 5-pt at 0,1
    d3[0, :] = (-5.0*wk[0, :] + 18.0*wk[1, :] - 24.0*wk[2, :] + 14.0*wk[3, :] - 3.0*wk[4, :]) / h3
    d3[1, :] = (-5.0*wk[1, :] + 18.0*wk[2, :] - 24.0*wk[3, :] + 14.0*wk[4, :] - 3.0*wk[5, :]) / h3
    # backward 5-pt at Nz-1, Nz-2
    d3[-1, :] = ( 5.0*wk[-1, :] - 18.0*wk[-2, :] + 24.0*wk[-3, :] - 14.0*wk[-4, :] + 3.0*wk[-5, :]) / h3
    d3[-2, :] = ( 5.0*wk[-2, :] - 18.0*wk[-3, :] + 24.0*wk[-4, :] - 14.0*wk[-5, :] + 3.0*wk[-6, :]) / h3
    return d3

# --------- bonus: tiny analytic check (optional) ---------
def _test_convergence():
    """
    Quick self-test: use w(z,x)=sin(m z)*cos(kx x) and check errors scale ~ O(dz^2).
    """
    m = 2.3
    Lx, Nx = 2*np.pi, 256
    Lz, Nz = 1.0, 128
    x = np.linspace(0, Lx, Nx, endpoint=False)
    z = np.linspace(0, Lz, Nz)
    dz = z[1]-z[0]
    X, Z = np.meshgrid(x, z)
    w = np.sin(m*Z)*np.cos(2*np.pi*x/Lx)

    # RFFT in x
    wk = np.fft.rfft(w, axis=1)

    d1_true =  m*np.cos(m*Z)*np.cos(2*np.pi*x/Lx)
    d2_true = -m*m*np.sin(m*Z)*np.cos(2*np.pi*x/Lx)
    d3_true = -m**3*np.cos(m*Z)*np.cos(2*np.pi*x/Lx)

    d1 = irfft(diff1_z_k(wk, dz), n=Nx, axis=1).real
    d2 = irfft(diff2_z_k(wk, dz), n=Nx, axis=1).real
    d3 = irfft(diff3_z_k(wk, dz), n=Nx, axis=1).real

    e1 = np.max(np.abs(d1 - d1_true))
    e2 = np.max(np.abs(d2 - d2_true))
    e3 = np.max(np.abs(d3 - d3_true))
    print(f"max errors: d1={e1:.2e}, d2={e2:.2e}, d3={e3:.2e}")
