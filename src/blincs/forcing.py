#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 09:39:26 2025

@author: Hossein Amini Kafiabad and Majid Bastankhah
(proudly collaborated with chatGPT to tidy up this module)
"""

# forcing.py
from __future__ import annotations
from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np
from scipy.fft import rfft
from scipy.integrate import cumulative_trapezoid


# ---------- data classes for the geometry of turbines/windfarms ----------
@dataclass
class BoxFarmGeom:
    ax: float
    az: float
    zh: float
    x0: float = 0.0   # center in x

    def as_rect(self):
        """Return rectangle bounds for plotting [xmin, xmax, zmin, zmax]."""
        return (self.x0 - self.ax/2, self.x0 + self.ax/2,
                self.zh - self.az/2, self.zh + self.az/2)



# ---------- utilities ----------
def _dfk_dz_from_Fxr(Fxr: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    Compute ∂f/∂z in real space, then RFFT in x → dfk_dz (Nz, Nk).
    """
    dz = float(z[1] - z[0])
    df_dz = np.gradient(Fxr, dz, axis=0)
    return rfft(df_dz, axis=1, workers=-1)


# ---------- base interface ----------
class ForcingBase(ABC):
    """
    Abstract base for forcing models.

    Contract: return (Fxr, dfk_dz) with shapes:
      Fxr:    (Nz, Nx)  real-space forcing field
      dfk_dz: (Nz, Nk)  RFFT in x of ∂Fxr/∂z
    """
    @abstractmethod
    def generate(
        self,
        U0z: np.ndarray,      # (Nz,)
        dU0z: np.ndarray,     # (Nz,)
        k_vec: np.ndarray,    # (Nk,)
        xx: np.ndarray,       # (Nz, Nx)
        zz: np.ndarray,       # (Nz, Nx)
        z:  np.ndarray,       # (Nz,)
        state: dict | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        ...


# ---------- concrete forcings ----------
@dataclass
class TurbineBoxForcing(ForcingBase):
    geom: BoxFarmGeom | None = None     # single geometry (optional)
    geoms: list[BoxFarmGeom] | None = None  # multiple farms (new feature)

    Fx: float | None = None
    D: float | None = None
    Ct: float | None = None
    Sx: float | None = None
    Sy: float | None = None
    smooth_x_factor: float = 100.0
    smooth_z_factor: float = 20.0

    def _single_forcing(self, geom, U0z, xx, zz, z):
        """Helper for one farm geometry."""
        ax, az, zh, x0 = geom.ax, geom.az, geom.zh, geom.x0
        
        
        # hub-height wind
        Uh = np.interp(zh, z, U0z)

        # forcing amplitude
        Fx = self.Fx if self.Fx is not None else (
            # -0.5 * self.Ct * Uh**2 * (np.pi * 0.25 * self.D**2) /
            # (self.Sx * self.Sy * az)
            -0.7 * self.Ct * Uh**2 * (np.pi * self.D) /
            (8*self.Sx * self.Sy)
        )
        
        # smoothed box forcing
        smooth_x = ax / self.smooth_x_factor
        smooth_z = az / self.smooth_z_factor
        return Fx / 4.0 * (
            (np.tanh((xx - x0)/smooth_x) - np.tanh((xx - x0 - ax)/smooth_x)) *
            (np.tanh((zz - zh + az/2)/smooth_z) - np.tanh((zz - zh - az/2)/smooth_z))
        )
    
  

    def generate(self, U0z, dU0z, k_vec, xx, zz, z, state=None):
        """Generate turbine forcing field in x–z space and its spectral derivative."""

        if self.geoms is not None:  # multiple farms
            Fxr = np.zeros_like(xx, dtype=float)
            for g in self.geoms:
                Fxr += self._single_forcing(g, U0z, xx, zz, z)
        elif self.geom is not None:  # single farm
            Fxr = self._single_forcing(self.geom, U0z, xx, zz, z)
        else:
            raise ValueError("Must provide either geom or geoms.")

        # spectral derivative in x
        dz_val = z[1] - z[0]
        df_dz = np.gradient(Fxr, dz_val, axis=0)
        dfk_dz = rfft(df_dz, axis=1, workers=-1)

        return Fxr, dfk_dz    
    



@dataclass
class GaussianForcing(ForcingBase):
    """Smooth Gaussian forcing in x–z."""
    amplitude: float
    x0: float
    z0: float
    sig_x: float
    sig_z: float

    def generate(self, U0z, dU0z, k_vec, xx, zz, z, state=None):
        Fxr = self.amplitude * np.exp(-((xx - self.x0)**2)/(2*self.sig_x**2)) \
                              * np.exp(-((zz - self.z0)**2)/(2*self.sig_z**2))
        dfk_dz = _dfk_dz_from_Fxr(Fxr, z)
        return Fxr, dfk_dz


@dataclass
class UserFieldForcing(ForcingBase):
    """
    Use a supplied real-space field Fxr_fn(x,z) or a prebuilt array.
    If array is provided, its shape must be (Nz, Nx).
    """
    Fxr_array: np.ndarray | None = None
    Fxr_fn: callable | None = None   # signature F(x_grid, z_grid) -> array (Nz, Nx)

    def generate(self, U0z, dU0z, k_vec, xx, zz, z, state=None):
        if self.Fxr_array is not None:
            Fxr = np.asarray(self.Fxr_array)
        elif self.Fxr_fn is not None:
            Fxr = np.asarray(self.Fxr_fn(xx, zz))
        else:
            raise ValueError("UserFieldForcing: provide either Fxr_array or Fxr_fn.")
        dfk_dz = _dfk_dz_from_Fxr(Fxr, z)
        return Fxr, dfk_dz







@dataclass
class NonlinearForcing(ForcingBase):
    """
    Nonlinear turbine forcing (single turbine), depending on the current flow.

    df/dz ∝ (U0 + ur)^2 * G_x(x; x0, sigma_x) * V_z(z; zh, R)

    where:
      - G_x is a Gaussian in x (no normalization, peak=1 at x=x0),
      - V_z is the vertical profile with edge regularization,
      - R = D/2 is turbine radius.

    Parameters
    ----------
    x0 : float
        Turbine center in x (m). Default 0.0.
    zh : float
        Hub height (m).
    D  : float
        Rotor/Turbine diameter (m) → R = D/2.
    Ct : float
        Ct (default 1.33).
    alpha : float
        multiplying factor (default 1.0).
    Sy : float | None
        Lateral spacing (m). If None → defaults to 4*D (not implemented yet).
    sigma_x : float | None
        Gaussian width in x (m). If None → defaults to D/2.

    reg_method : {"cos_taper", "bump", "no_reg"}
        Choice of regularisation for vertical profile
        
    beta : float
        For "smooth bump": Sharpness at |z-zh|=R (larger -> more damped on the edges).
    eps_offset : float
        For cosine tapering: the offset from the edges that is smoothed out
    """ 
   

    x0: float = 0.0
    zh: float = 100.0
    D: float = 100.0
    Ct: float = 1.33
    alpha: float = 1.0
    Sy: float | None = None
    sigma_x: float | None = None
    reg_method: str = "bump"
    beta: float = 0.03
    eps_offset: float = 0.05

    # --- vertical profile with regularised edges ---
    def _edge_shapes(self, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Rotor-thickness profile derivative, regularised near |s|=R.
        Returns V(z) and S(z) of shape (Nz,).
        """
        R = 0.5 * self.D
        s = (z - self.zh) / R
        S = np.zeros_like(s, dtype=float)
        V = np.zeros_like(s, dtype=float)
        inside = np.abs(s) < 1.0
        if not np.any(inside):
            return S, V

        s_in = s[inside]
        base_S = R * np.sqrt(1.0 - s_in**2)         # sqrt(R^2 - (z-zh)^2)
        base_V = (R * s_in) / base_S                # (z-zh)/sqrt(...)

        if self.reg_method == "cos_taper":
            taper = np.ones_like(s_in)
            band = (np.abs(s_in) > (1.0 - self.eps_offset))
            if np.any(band):
                rel = (np.abs(s_in[band]) - (1.0 - self.eps_offset)) / self.eps_offset
                taper[band] = 0.5 * (1.0 + np.cos(np.pi * rel))
            base_S *= taper
            base_V *= taper
        elif self.reg_method == "bump":
            bump = np.exp(-self.beta * s_in**2 / (1.0 - s_in**2))
            base_S *= bump
            base_V *= bump
        elif self.reg_method == "no_reg":
            pass
        else:
            raise ValueError("Unknown reg_method; use 'cos_taper' or 'bump'.")

        S[inside] = base_S
        V[inside] = base_V
        return S, V
    
    @staticmethod
    def _dz_centered(arr: np.ndarray, z: np.ndarray) -> np.ndarray:
        Nz, Nx = arr.shape
        dz = z[1] - z[0]
        out = np.empty_like(arr)
        out[1:-1, :] = (arr[2:, :] - arr[:-2, :]) / (2.0 * dz)
        out[0,  :] = (-3*arr[0, :] + 4*arr[1, :] - arr[2, :]) / (2.0 * dz)
        out[-1, :] = ( 3*arr[-1, :] - 4*arr[-2, :] + arr[-3, :]) / (2.0 * dz)
        return out
  
    # --- main forcing ---
    def generate(self, U0z, dU0z, k_vec, xx, zz, z, state=None):
        """
        Build df/dz and F_x with corrected signs:
            df/dz = (Ct α^2 / Sy) * Gx(x) * [ (U0+u)^2 V(z) - 2(U0'+u_z)(U0+u) S(z) ].

        Parameters
        ----------
        U0z : (Nz,) background profile
        dU0z: (Nz,) background vertical derivative (explicit input)
        state: dict with 'ur': (Nz,Nx) current u-perturbation

        Returns
        -------
        Fxr   : (Nz,Nx) real
        dfk_dz: (Nz,Nk) complex (RFFT in x of df/dz)
        """
        if state is None or ("ur" not in state):
            ur = np.zeros_like(xx)
        else:
            ur = np.asarray(state["ur"])         # (Nz,Nx)

        Sy = float(4.0 * self.D) if (self.Sy is None) else float(self.Sy)
        sigma_x = (self.D / 2.0) if (self.sigma_x is None) else float(self.sigma_x)

        # x-localization (unnormalized Gaussian, peak=1 at x0)
        x_row = xx[0, :]
        Gx = np.exp(-0.5 * ((x_row - self.x0) / sigma_x) ** 2)

        # regularized edge shapes
        S_z, V_z = self._edge_shapes(z)

        # total wind and vertical derivative
        Utot = U0z[:, None] + ur            # (Nz,Nx)
        du_dz = self._dz_centered(ur, z)    # (Nz,Nx)
        dUtot_z = dU0z[:, None] + du_dz      # (Nz,Nx)

        # amplitude
        amp = (self.Ct * (self.alpha ** 2)) / Sy

        # corrected signs:
        df_dz = amp * (
            (Utot**2) * V_z[:, None]
            - 2.0 * dUtot_z * Utot * S_z[:, None]
        ) * Gx[None, :]

        # integrate in z → F_x (zero at bottom)
        Fxr = cumulative_trapezoid(df_dz, z, axis=0, initial=0.0)

        # to spectral-x
        dfk_dz = rfft(df_dz, axis=1, workers=-1)

        return Fxr, dfk_dz





# ---------- simple factory ----------
def make_forcing(kind: str, **kwargs) -> ForcingBase:
    kind = kind.lower()
    if kind in ("turbine", "turbine_box", "box"):
        return TurbineBoxForcing(**kwargs)
    if kind in ("gaussian", "gauss"):
        return GaussianForcing(**kwargs)
    if kind in ("user", "user_field", "custom"):
        return UserFieldForcing(**kwargs)
    if kind in ("nonlinear", "feedback"):
        return NonlinearForcing(**kwargs)
    raise ValueError(f"Unknown forcing kind: {kind}")

