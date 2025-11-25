# Example 1 — Wave‑forced Boundary Layer

```python
from blincs.blwave_simulation import simulate_BLwave
from blincs.forcing import TurbineBoxForcing, BoxFarmGeom
from blincs.vel_base import LogProfile

forcing = TurbineBoxForcing(BoxFarmGeom(ax=20000, az=100, zh=100))
background = LogProfile(U_top=20.0, z0=1.0, ztop=1000.0)

x, z, ur, wr, eta, Fxr, U0z = simulate_BLwave(
    g_prime=0.3, N=0.01,
    background=background, forcing=forcing,
    Nx=1024, Nz=200, Lx=8_000_000, Lz=1000
)
```
