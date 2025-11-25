#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
The main script where the simulation are run and plotted

'''
import numpy as np

# --- import both versions ---
from BLwave_simulation import simulate_BLwave
from U_background import LogProfile
from farm_turbine_forcing import ForcingBase, NonlinearForcing, BoxFarmGeom, TurbineBoxForcing
from plotting import plot_fields_fixedz, compare_fields_fixedz, plot_field_xz
import time

# ---------------- Parameters ----------------
Nx = 1024*8
Nz = 1000
Lx = 8000000
Lz = 1000
Uwind = 20.0
g_prime = 0.3
N = 0.01
nu_turb = 10.0

# geometry
farm_geom = BoxFarmGeom(ax=20000, az=100, zh=100, x0=0.0)

# forcing setup
forcing = TurbineBoxForcing(
    geom=farm_geom,
    Fx=None,   # let turbine physics compute it
    D=100.0,
    Ct=1.5,
    Sx=1000.0,
    Sy=1000.0,
)

# --- forcing instance ---
# forcing = NonlinearForcing(
#     x0=0.0,
#     zh=100.0,
#     D=100.0,
#     Ct=.02,
#     alpha=1.0,
#     Sy=400.,
#     sigma_x=100.0,       # narrow Gaussian in x
#     reg_method="bump",  # try also "bump"
#     beta=0.1,
#     eps_offset=0.1,
# )


# background profile
background = LogProfile(U_top=Uwind, z0=1.0, ztop=Lz)

t0 = time.time()
x, z, ur, wr, eta, Fxr, U0z = simulate_BLwave(
    g_prime=g_prime, N=N,
    background=background, forcing=forcing,
    Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, Uwind=Uwind, nu_turb=nu_turb,
    BC_method="one-sided",
    nonlinear_adv = True,
    max_iter = 60,
    tol = 1e-5,
    dealias_NL = True,
)
print('simulation time = ' , time.time()-t0)

########## PLOTS ################

fields = {
    "u_r": ur,
    "w_r": wr,
    "eta": eta,   # 1D → automatically handled
}

# Single z-target
plot_fields_fixedz(
    x, z, fields, z_target=100,
    geom=farm_geom,
    units={"u_r":"m/s", "w_r":"m/s", "eta":"m"},
    xlim=(-50,200)
)

# Multiple z-targets
plot_fields_fixedz(
    x, z, {"u_r": ur}, z_target=[50, 100, 500], xlim=(-100, 200),
    geom=farm_geom,
    units={"u_r":"m/s"}
)

# # ur as contourf, symmetric limits
# plot_field_xz(x, z, ur, field_name="u_r", units="m/s",
#               mode="contourf", levels=40, symmetric=True,
#               farm_geom=farm_geom, xlim=(-50, 200))

# # wr as pcolormesh for debugging (no interpolation)
# plot_field_xz(x, z, wr, field_name="w_r", units="m/s", xlim=(-50, 200),
#               mode="pcolormesh", cmap="RdBu_r",
#               farm_geom=farm_geom)

x2, z2, ur2, wr2, eta2, Fxr2, U0z2 = simulate_BLwave(
    g_prime=g_prime, N=N,
    background=background, forcing=forcing,
    Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, Uwind=Uwind, nu_turb=nu_turb,
    BC_method="one-sided derivative",
    nonlinear_adv = False
)

res1 = {"ur": ur, "wr": wr, "eta": eta, "z": z, "x": x}
res2 = {"ur": ur2, "wr": wr2, "eta": eta2, "z": z2, "x": x2}

# # compare ur and eta together
compare_fields_fixedz(
    [res1, res2],
    fields=["ur", "eta"],
    z_target=100.0,
    labels=["nonlinear", "linear"],
    farm_geom=farm_geom,
    units={"ur": "m/s", "eta": "m"},
    xlimits=(-50, 200)
)

compare_fields_fixedz(
    [res1, res2],
    fields=["ur", "wr"],
    z_target=400.0,
    labels=["nonlinear", "linear"],
    farm_geom=farm_geom,
    units={"ur": "m/s", "wr": "m/s"},
    xlimits=(-50, 200)
)




# # comparing two or more separate simulations with one another
# x2, z2, ur2, wr2, eta2, Fxr2, U0z2 = simulate_BLwave(
#     g_prime=g_prime, N=N,
#     background=background, forcing=forcing,
#     Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, Uwind=Uwind, nu_turb=nu_turb,
#     BC_method="one-sided derivative"
# )

# Nx = 1024*10
# Nz = 1400

# x3, z3, ur3, wr3, eta3, Fxr3, U0z3 = simulate_BLwave(
#     g_prime=g_prime, N=N,
#     background=background, forcing=forcing,
#     Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, Uwind=Uwind, nu_turb=nu_turb,
#     BC_method="ghost point"
# )

# res1 = {"ur": ur, "wr": wr, "eta": eta, "z": z, "x": x}
# res2 = {"ur": ur2, "wr": wr2, "eta": eta2, "z": z2, "x": x2}
# res3 = {"ur": ur3, "wr": wr3, "eta": eta3, "z": z3, "x": x3}

# # compare ur and eta together
# compare_fields_fixedz(
#     [res1, res2, res3],
#     fields=["ur", "eta"],
#     z_target=100.0,
#     labels=["ghost-point", "one-sided derivative", "high resolution"],
#     farm_geom=farm_geom,
#     units={"ur": "m/s", "eta": "m"},
#     xlimits=(-50, 200)
# )





