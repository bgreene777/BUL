'''
Script to produce quick time-height figures of data from AERIoe retrievals
using data from the DOE ARM SGP site in spring 2017
Dates analyzed: 20--23 May

Author: Brian Greene
University of Oklahoma
Last update: February 2020
'''
# Python Packages
import os
from datetime import datetime, timedelta
import warnings

# Installed packages
import netCDF4
import numpy as np
import cmocean

from glob import glob
from matplotlib import rc
from matplotlib import dates as mpdates
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# ---------------- #
# Data directories #
# ---------------- #
# home directory
home = os.path.expanduser("~")
# main data directory
BUL = os.path.join(home, "Documents", "Data", "BUL")
# figure save directory
figpath = os.path.join(BUL, "Figures", "Quicklook")
if not os.path.exists(figpath):
    os.mkdir(figpath)

# -------------------------- #
# Load Water Vapor Dial Data #
# -------------------------- #
# AERIonly
f_AERI = glob(os.path.join(BUL, "WaterVaporDial", "AERIonly", "*.2017052[0-3]*.cdf"))
AERI_dic = {}
d = 20
# loop through these files and load all vars and ncattrs into dictionary
print("---AERI Only---")
for f in f_AERI:
    AERI_dic[d] = {}
    print(f"Reading file: {f.split(os.sep)[-1]}")
    df = netCDF4.Dataset(f, "r")
    for var in df.variables:
        AERI_dic[d][var] = df.variables[var][:]
    for attr in df.ncattrs():
        AERI_dic[d][attr] = df.getncattr(attr)
    df.close()
    d += 1

#AERIrLID
print("---AERI with Raman Lidar---")
f_AERIrLID = glob(os.path.join(BUL, "WaterVaporDial", "AERIrLID", "*.2017052[0-3]*.cdf"))
AERIrLID_dic = {}
d = 20
# loop through these files and load all vars and ncattrs into dictionary
for f in f_AERIrLID:
    AERIrLID_dic[d] = {}
    print(f"Reading file: {f.split(os.sep)[-1]}")
    df = netCDF4.Dataset(f, "r")
    for var in df.variables:
        AERIrLID_dic[d][var] = df.variables[var][:]
    for attr in df.ncattrs():
        AERIrLID_dic[d][attr] = df.getncattr(attr)
    df.close()
    d += 1

#AERIvDIAL
print("---AERI with Water Vapor Dial---")
f_AERIvDIAL = glob(os.path.join(BUL, "WaterVaporDial", "AERIvDIAL", "*.2017052[0-3]*.cdf"))
AERIvDIAL_dic = {}
d = 20
# loop through these files and load all vars and ncattrs into dictionary
for f in f_AERIvDIAL:
    AERIvDIAL_dic[d] = {}
    print(f"Reading file: {f.split(os.sep)[-1]}")
    df = netCDF4.Dataset(f, "r")
    for var in df.variables:
        AERIvDIAL_dic[d][var] = df.variables[var][:]
    for attr in df.ncattrs():
        AERIvDIAL_dic[d][attr] = df.getncattr(attr)
    df.close()
    d += 1

# ------------------------------- #
# Plot Time-Height Cross Sections #
# ------------------------------- #
rc('font',weight='normal',size=20,family='serif',serif='Computer Modern Roman')
rc('text',usetex='True')
# AERI only
# create meshgrid
for d in range(20, 24):
    iz = np.where(AERI_dic[d]["height"] <= 4.)[0]
    X1, Y1 = np.meshgrid(AERI_dic[d]["hour"], AERI_dic[d]["height"][iz])

    fig1, ax1 = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(16,12))

    cfax1 = ax1[0].pcolormesh(X1, Y1, AERI_dic[d]["temperature"][:, iz].transpose(), 
        cmap=cmocean.cm.thermal, vmin=0., vmax=30.)
    cfax1.set_edgecolor("face")
    cbar1 = plt.colorbar(cfax1, ax=ax1[0])
    cbar1.ax.set_ylabel("Temperature [$^\circ$C]")
    # ax1[0].set_xlabel("Hour [UTC]")
    ax1[0].set_ylabel("Altitude [km AGL]")
    # ax1[0].set_xlim([0, 24])
    ax1[0].set_ylim([0, 4])
    ax1[0].yaxis.set_major_locator(MultipleLocator(1))
    ax1[0].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax1[0].set_title(f"2017-05-{d} AERI Only")

    # AERI + Raman
    iz = np.where(AERIrLID_dic[d]["height"] <= 4.)[0]
    X2, Y2 = np.meshgrid(AERIrLID_dic[d]["hour"], AERIrLID_dic[d]["height"][iz])

    cfax2 = ax1[1].pcolormesh(X2, Y2, AERIrLID_dic[d]["temperature"][:, iz].transpose(), 
        cmap=cmocean.cm.thermal, vmin=0., vmax=30.)
    cfax2.set_edgecolor("face")
    cbar2 = plt.colorbar(cfax2, ax=ax1[1])
    cbar2.ax.set_ylabel("Temperature [$^\circ$C]")
    # ax1.set_xlabel("Hour [UTC]")
    ax1[1].set_ylabel("Altitude [km AGL]")
    # ax1.set_xlim([0, 24])
    ax1[1].set_ylim([0, 4])
    ax1[1].yaxis.set_major_locator(MultipleLocator(1))
    ax1[1].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax1[1].set_title("AERI + Raman")

    # AERI + vDial
    iz = np.where(AERIvDIAL_dic[d]["height"] <= 4.)[0]
    X3, Y3 = np.meshgrid(AERIvDIAL_dic[d]["hour"], AERIvDIAL_dic[d]["height"][iz])

    cfax3 = ax1[2].pcolormesh(X3, Y3, AERIvDIAL_dic[d]["temperature"][:, iz].transpose(), 
        cmap=cmocean.cm.thermal, vmin=0., vmax=30.)
    cfax3.set_edgecolor("face")
    cbar3 = plt.colorbar(cfax3, ax=ax1[2])
    cbar3.ax.set_ylabel("Temperature [$^\circ$C]")
    ax1[2].set_xlabel("Hour [UTC]")
    ax1[2].set_ylabel("Altitude [km AGL]")
    ax1[2].set_xlim([0, 24])
    ax1[2].set_ylim([0, 4])
    ax1[2].xaxis.set_major_locator(MultipleLocator(3))
    ax1[2].xaxis.set_minor_locator(MultipleLocator(1))
    ax1[2].yaxis.set_major_locator(MultipleLocator(1))
    ax1[2].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax1[2].set_title("AERI + vDIAL")

    fig1.tight_layout()

    fig1.savefig(os.path.join(figpath, f"201705{d}_Temperature.pdf"), dpi=150, fmt="pdf")
    plt.close(fig1)

# ---------------------------------------------------------- #
# Plot differences of each time-height relative to AERI Only #
# ---------------------------------------------------------- #
# grab data
z = AERI_dic[20]["height"]
hrs = AERI_dic[20]["hour"]
iz = np.where(z <= 4.)[0]
it = np.where(AERIrLID_dic[20]["hour"] in hrs)[0]
AERI_T = AERI_dic[20]["temperature"][:, iz].transpose()
AERIrLID_T = AERIrLID_dic[20]["temperature"][:, iz][it].transpose()
AERIvDIAL_T = AERIvDIAL_dic[20]["temperature"][:, iz][it].transpose()
# calculate differences
AERI_LID_T_diff = AERI_T - AERIrLID_T
AERI_DIAL_T_diff = AERI_T - AERIvDIAL_T
# plot
fig2, ax2 = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(16, 8))
# AERI - AERIrLID
cfax21 = ax2[0].pcolormesh(hrs, z, AERI_LID_T_diff)

cfax22 = ax2[1].pcolormesh(hrs, z, AERI_DIAL_T_diff)