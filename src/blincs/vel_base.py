#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# U_background.py
from __future__ import annotations
from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np


class BackgroundBase(ABC):
    """Return (U0, dU0, d2U0), each shape (Nz,)."""
    @abstractmethod
    def eval(self, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...


# ---------------- Profiles ----------------

@dataclass
class ConstantFlow(BackgroundBase):
    """U0(z) = U (constant)."""
    U: float
    def eval(self, z: np.ndarray):
        z = np.asarray(z, dtype=float)
        Nz = z.size
        U0   = np.full(Nz, self.U, dtype=float)
        dU0  = np.zeros(Nz, dtype=float)
        d2U0 = np.zeros(Nz, dtype=float)
        return U0, dU0, d2U0


@dataclass
class LogProfile(BackgroundBase):
    """
    Logarithmic profile anchored by U_top at z=ztop:
        U0(z) = U_top * log(z / z0) / log(ztop / z0)
    Values at z <= z0 are clamped slightly above z0 to avoid singularities.
    """
    U_top: float
    z0: float = 1.0
    ztop: float | None = None

    def eval(self, z: np.ndarray):
        z = np.asarray(z, dtype=float)
        ztop = float(z[-1]) if self.ztop is None else float(self.ztop)
        zc   = np.maximum(z, self.z0 * (1.0 + 1e-12))
        denom = np.log(ztop / self.z0)
        U0   = self.U_top * np.log(zc / self.z0) / denom
        dU0  = self.U_top / (denom * zc)
        d2U0 = -self.U_top / (denom * zc**2)
        return U0, dU0, d2U0


@dataclass
class ParabolicTopCap(BackgroundBase):
    """
    Quadratic profile with a flat top (laminar-BL surrogate):

      Constraints:
        U0(z_min)   = U_bottom
        U0(z_max)   = U_top
        dU0/dz|z_max = 0

      Let z_min = z[0], z_top = z[-1]. Solve U = a z^2 + b z + c with:
        b = -2 a z_top
        a = (U_bottom - U_top) / (z_min - z_top)^2
        c = U_top + a z_top^2
    """
    U_top: float
    U_bottom: float

    def eval(self, z: np.ndarray):
        z = np.asarray(z, dtype=float)
        if not np.all(np.diff(z) > 0):
            raise ValueError("ParabolicTopCap: z must be strictly increasing.")
        z_min = float(z[0])
        z_top = float(z[-1])
        denom = (z_min - z_top)
        if abs(denom) < 1e-14:
            raise ValueError("ParabolicTopCap: z_min and z_top coincide.")
        a = (self.U_bottom - self.U_top) / (denom * denom)
        b = -2.0 * a * z_top
        c = self.U_top + a * (z_top ** 2)

        U0   = a * z**2 + b * z + c
        dU0  = 2.0 * a * z + b
        d2U0 = np.full_like(z, 2.0 * a, dtype=float)
        return U0, dU0, d2U0


@dataclass
class UserProfile(BackgroundBase):
    """User-supplied callable: fn(z) -> (U0, dU0, d2U0)."""
    fn: callable
    def eval(self, z: np.ndarray):
        return self.fn(np.asarray(z, dtype=float))


# --------------- Factory ----------------

def make_background(kind: str, **kwargs) -> BackgroundBase:
    k = kind.lower()
    if k in ("const", "constant"):
        return ConstantFlow(**kwargs)
    if k in ("log", "logarithmic"):
        return LogProfile(**kwargs)
    if k in ("parabolic", "parabola", "parabolic_top", "parab_top"):
        return ParabolicTopCap(**kwargs)
    if k in ("user", "custom"):
        return UserProfile(**kwargs)
    raise ValueError(f"Unknown type of background velocity: {kind}")
