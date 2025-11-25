#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLOTTTING FUNCTIONS
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from typing import Sequence, List, Tuple, Optional, Union, Literal
from .forcing import BoxFarmGeom



def plot_fields_fixedz(
    x: np.ndarray,
    z: np.ndarray,
    fields: dict[str, np.ndarray],
    z_target: float | Sequence[float],
    geom: BoxFarmGeom | None = None,
    xlim: tuple[float, float] | None = None,
    normalize_forcing: tuple[np.ndarray, str] | None = None,
    units: dict[str, str] | None = None,
    figsize=(8, 6)
):
    """
    Plot one or more fields at one or more fixed heights.

    Parameters
    ----------
    x : array (Nx,)
        Horizontal coordinate.
    z : array (Nz,)
        Vertical coordinate.
    fields : dict[str, array]
        {name: field}, each of shape (Nz, Nx) or (Nx,) for 1D fields like eta.
    z_target : float or list of floats
        Height(s) at which to extract slices. Ignored for 1D fields.
    geom : BoxFarmGeom, optional
        Farm geometry for vertical markers.
    xlim : (xmin, xmax), optional
        X-axis limits.
    normalize_forcing : (Fxr, label), optional
        If given, overlays scaled forcing profile at z_target (only if z_target is scalar).
    units : dict[str, str], optional
        Mapping {field_name: units}.
    figsize : tuple
        Figure size.

    Returns
    -------
    fig, axes
    """
    # Ensure z_target is iterable
    if np.isscalar(z_target):
        z_targets = [float(z_target)]
    else:
        z_targets = list(z_target)

    nplots = len(fields)
    fig, axes = plt.subplots(nplots, 1, figsize=figsize, sharex=True)
    if nplots == 1:
        axes = [axes]

    for ax, (name, field) in zip(axes, fields.items()):
        field = np.asarray(field)

        if field.ndim == 2:  # (Nz, Nx) field
            for zt in z_targets:
                iz = int(np.argmin(np.abs(z - zt)))
                ax.plot(x/1000, field[iz, :], label=f"{name}(z={zt:.1f}m)")
        elif field.ndim == 1:  # (Nx,) field (like eta)
            ax.plot(x/1000, field, label=f"{name}")
        else:
            raise ValueError(f"Field {name} has unsupported shape {field.shape}")

        ax.set_ylabel(f"{name} [{units.get(name,'') if units else ''}]")
        ax.grid(True)

        # overlay forcing if given AND if z_target is scalar
        if normalize_forcing is not None and len(z_targets) == 1 and field.ndim == 2:
            Fxr, flabel = normalize_forcing
            iz = int(np.argmin(np.abs(z - z_targets[0])))
            overlay = -0.5 * np.abs(Fxr[iz, :]) / np.max(np.abs(Fxr[iz, :]))
            ax.plot(x/1000.0, overlay, 'k--', lw=1, label=flabel)

        # vertical lines for farm geometry
        if geom is not None:
            ax.axvline((geom.x0 - geom.ax/2) / 1000, color='r', ls='--', lw=0.8)
            ax.axvline((geom.x0 + geom.ax/2) / 1000, color='r', ls='--', lw=0.8)

        if xlim is not None:
            ax.set_xlim(xlim)

        ax.legend()

    axes[-1].set_xlabel("x [km]")
    fig.tight_layout()
    return fig, axes



def compare_fields_fixedz(
    results: List[dict],
    fields: Union[str, List[str]],
    z_target: float,
    ylabel_name: Optional[str] = None,
    labels: Optional[List[str]] = None,
    farm_geom: Optional[object] = None,  # BoxFarmGeom or similar
    xlimits: Optional[Tuple[float, float]] = None,
    units: Optional[dict] = None,
    km: bool = True,
    figsize: Tuple[float, float] = (8, 3),
):
    """
    Compare one or more fields at a fixed height across multiple simulations,
    each with its own horizontal resolution.

    Parameters
    ----------
    results : list of dict
        Each dict must contain keys "x" and the requested fields.
        If a field is 2D (Nz,Nx), the dict must also include "z".
    fields : str or list of str
        Field(s) to compare (e.g. "ur", "wr", "eta").
    z_target : float
        Height (m) to extract. Ignored if field is 1D (Nx,).
    labels : list of str, optional
        Legend labels for each result. Defaults to "Run 0", "Run 1", ...
    farm_geom : BoxFarmGeom, optional
        If provided, vertical dashed lines mark farm edges.
    xlimits : (xmin, xmax), optional
        Limits for x-axis; same units as plotted x (km if km=True else m).
    units : dict, optional
        Mapping field → units string for y-axis labels.
        Example: {"ur": "m/s", "eta": "m"}.
    km : bool
        If True, plot x in km; otherwise in meters.
    figsize : tuple
        Figure size *per subplot*. Height is scaled automatically.

    Returns
    -------
    fig, axes : matplotlib Figure and list of Axes
    """
    if isinstance(fields, str):
        fields = [fields]

    if labels is None:
        labels = [f"Run {i}" for i in range(len(results))]

    fig, axes = plt.subplots(
        len(fields), 1,
        figsize=(figsize[0], figsize[1]*len(fields)),
        sharex=False
    )
    if len(fields) == 1:
        axes = [axes]

    for ax, field in zip(axes, fields):
        for res, lab in zip(results, labels):
            x = np.asarray(res["x"])
            x_plot = x / 1000.0 if km else x
            arr = np.asarray(res[field])

            if arr.ndim == 2:
                if "z" not in res:
                    raise ValueError(f'result dict must include "z" for 2D field "{field}".')
                z = np.asarray(res["z"])
                iz = int(np.argmin(np.abs(z - z_target)))
                ax.plot(x_plot, arr[iz, :], label=f"{lab}")
            elif arr.ndim == 1:
                ax.plot(x_plot, arr, label=lab)
            else:
                raise ValueError(f'Field "{field}" has unsupported shape {arr.shape}.')

        # farm edges (draw once, all runs assumed same farm_geom)
        if farm_geom is not None:
            x1 = (farm_geom.x0 - farm_geom.ax/2) / (1000.0 if km else 1.0)
            x2 = (farm_geom.x0 + farm_geom.ax/2) / (1000.0 if km else 1.0)
            ax.axvline(x1, color="r", ls="--", lw=0.9)
            ax.axvline(x2, color="r", ls="--", lw=0.9)

        if xlimits is not None:
            ax.set_xlim(xlimits)

        ylabel = field
        if units and field in units:
            ylabel += f" [{units[field]}]"
        if ylabel_name is not None:
            ylabel = ylabel_name
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend()

    axes[-1].set_xlabel("x [km]" if km else "x [m]")
    fig.tight_layout()
    return fig, axes




def plot_field_xz(
    x: np.ndarray,
    z: np.ndarray,
    field: np.ndarray,                 # shape (Nz, Nx)
    field_name: str = "field",
    units: Optional[str] = None,
    *,
    mode: Literal["contourf", "pcolormesh"] = "contourf",
    levels: int | np.ndarray = 40,
    cmap: str = "RdBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    symmetric: bool = False,           # if True, enforce vmin=-vmax
    extend: Literal["neither","both","min","max"] = "neither",
    farm_geom: Optional[BoxFarmGeom] = None,
    xlim: Optional[Tuple[float, float]] = None,  # in km (if km=True) or m
    zlim: Optional[Tuple[float, float]] = None,
    km: bool = True,                   # plot x in km if True
    figsize: Tuple[float, float] = (10, 4),
    add_colorbar: bool = True,
):
    """
    Plot a 2D field in the (x,z) plane with either contourf or pcolormesh.

    Parameters
    ----------
    x, z : 1D arrays
        Grids in meters (Nx,), (Nz,).
    field : 2D array
        Data on (Nz, Nx).
    field_name : str
        For title / colorbar label.
    units : str, optional
        For colorbar label (e.g. 'm/s').
    mode : 'contourf' or 'pcolormesh'
        Rendering style.
    levels : int or array-like
        Number of levels (int) or explicit level values for contourf.
    cmap : str
        Matplotlib colormap name.
    vmin, vmax : floats, optional
        Color limits. If None, inferred from data.
    symmetric : bool
        If True, use v = max(|min|, |max|) and set (vmin, vmax)=(-v, v).
    extend : str
        Colorbar extend option for contour plots.
    farm_geom : BoxFarmGeom, optional
        If provided, draws the farm rectangle.
    xlim, zlim : (min, max), optional
        Axis limits (km if km=True else meters).
    km : bool
        Plot x in km if True; else in meters.
    """
    field = np.asarray(field)
    assert field.shape == (len(z), len(x)), "field must be (Nz, Nx)"

    # axes units
    x_plot = x / 1000.0 if km else x
    xlab = "x [km]" if km else "x [m]"

    # color limits
    if vmin is None or vmax is None:
        fmin, fmax = np.nanmin(field), np.nanmax(field)
        if symmetric:
            v = max(abs(fmin), abs(fmax))
            vmin = -v if vmin is None else vmin
            vmax =  v if vmax is None else vmax
        else:
            vmin = fmin if vmin is None else vmin
            vmax = fmax if vmax is None else vmax

    # mesh
    XX, ZZ = np.meshgrid(x_plot, z)

    fig, ax = plt.subplots(figsize=figsize)

    if mode == "contourf":
        cs = ax.contourf(XX, ZZ, field, levels=levels, cmap=cmap,
                         vmin=vmin, vmax=vmax, extend=extend)
    elif mode == "pcolormesh":
        # pcolormesh expects bin edges for perfect alignment; for simplicity use centers
        cs = ax.pcolormesh(XX, ZZ, field, cmap=cmap, shading="auto",
                           vmin=vmin, vmax=vmax)
    else:
        raise ValueError("mode must be 'contourf' or 'pcolormesh'")

    # farm rectangle
    if farm_geom is not None:
        x1 = (farm_geom.x0 - farm_geom.ax/2) / (1000.0 if km else 1.0)
        y1 = farm_geom.zh - farm_geom.az/2
        width = (farm_geom.ax) / (1000.0 if km else 1.0)
        height = farm_geom.az
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((x1, y1), width, height,
                               fill=False, ec="k", lw=1.2))

    if xlim is not None:
        ax.set_xlim(xlim)
    if zlim is not None:
        ax.set_ylim(zlim)

    ax.set_xlabel(xlab)
    ax.set_ylabel("z [m]")
        
    ax.grid(False)

    if add_colorbar:
        cbar = fig.colorbar(cs, ax=ax, pad=0.02)
        if units:
            cbar.set_label(units)
        tick_locator = ticker.MaxNLocator(nbins=4)
        cbar.locator = tick_locator
        cbar.update_ticks()

    fig.tight_layout()
    return fig, ax

