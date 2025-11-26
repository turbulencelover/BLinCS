#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script runs simulations for multiple farms using the BLINCS framework.
It sets up the simulation parameters, defines the farm geometry and forcing,
and executes the simulation to obtain the flow fields.

'''

import numpy as np
import matplotlib.pyplot as plt
import blincs
# --- import both versions ---
from blincs.blwave_simulation import simulate_BLwave
from blincs.vel_base import LogProfile
from blincs.forcing import ForcingBase, NonlinearForcing, BoxFarmGeom, TurbineBoxForcing
from blincs.plotting import plot_fields_fixedz, compare_fields_fixedz, plot_field_xz

# ---------------- Parameters ----------------
Nx = 1024*8
Nz = 1000
Lx = 2000000
Lz = 1000
Uwind = 20.0
g_prime = 0
N = 0.
nu_turb = 2.0


# ==== farm geometry and forcing =====
# we put two farms in the domain (separated by 9 farm lengths)
farm_geom = [
    BoxFarmGeom(ax=100.*105., az=100, zh=100, x0=0),
    BoxFarmGeom(ax=100.*105., az=100, zh=100, x0=+9*100.*105+100.*105.),
    ]

forcing = TurbineBoxForcing(
    geoms=farm_geom,
    Fx=None,   # let turbine physics compute it
    D=100.0,
    Ct=.8,
    Sx=7.*100.,
    Sy=5.*100.,
)

# ==== Background (Base) velocity =====
background = LogProfile(U_top=Uwind, z0=1.0, ztop=Lz)

# --- run simulation ---
x1, z1, u1, wr1, eta1, Fxr1, U0z1 = simulate_BLwave(
    g_prime=g_prime, N=N,
    background=background, forcing=forcing,
    Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, Uwind=Uwind, nu_turb=nu_turb,
    BC_method="one-sided",
    nonlinear_adv = False,
)

# --- plotting ---
# find reference velocity at hub height to normalize velocity deficit    
zh =100
z = np.linspace(1.0, Lz, Nz)
U0z, dU0z, d2U0z = background.eval(z)
Uh = np.interp(zh, z, U0z)

plot_field_xz(x1, z1, u1/Uh, field_name="u",
              mode="contourf", levels=40, symmetric=False,
              cmap = "bone", vmin= -.6, vmax = 0,
              figsize = (6, 1.5), xlim=(-50,450), zlim=(0,800)
              )

plot_fields_fixedz(x1, z1, {"u_r": u1/Uh}, 
                z_target=[50, 100, 500, 700], xlim=(-100, 400),
                units={"u_r":"m/s"})
        
plt.show()