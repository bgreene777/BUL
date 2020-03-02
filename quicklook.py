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
from scipy.interpolate import interp1d
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
# sort files
date_AERI = [int(ff.split(os.sep)[-1].split(".")[-3]) for ff in f_AERI]
id_AERI = np.argsort(date_AERI)
AERI_dic = {}
d = 20
# loop through these files and load all vars and ncattrs into dictionary
print("---AERI Only---")
for f in np.asarray(f_AERI)[id_AERI]:
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
date_AERIrLID = [int(ff.split(os.sep)[-1].split(".")[-3]) for ff in f_AERIrLID]
id_AERIrLID = np.argsort(date_AERIrLID)
AERIrLID_dic = {}
d = 20
# loop through these files and load all vars and ncattrs into dictionary
for f in np.asarray(f_AERIrLID)[id_AERIrLID]:
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
date_AERIvDIAL = [int(ff.split(os.sep)[-1].split(".")[-3]) for ff in f_AERIvDIAL]
id_AERIvDIAL = np.argsort(date_AERIvDIAL)
AERIvDIAL_dic = {}
d = 20
# loop through these files and load all vars and ncattrs into dictionary
for f in np.asarray(f_AERIvDIAL)[id_AERIvDIAL]:
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
print("---Plotting---")
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
    fsave1 = f"201705{d}_Temperature.pdf"
    print(f"Saving figure: {fsave1}")
    fig1.savefig(os.path.join(figpath, fsave1), dpi=150, fmt="pdf")
    plt.close(fig1)

# ---------------------------------------------------------- #
# Plot differences of each time-height relative to AERI Only #
# ---------------------------------------------------------- #
for d in range(20, 24):
    # grab data
    # heights below 4 km (all files symmetric)
    z = AERI_dic[d]["height"]
    iz = np.where(z <= 4.)[0]
    # missing temporal data in some files
    # find longest temporal array
    hrs1 = AERI_dic[d]["hour"]
    hrs2 = AERIrLID_dic[d]["hour"]
    hrs3 = AERIvDIAL_dic[d]["hour"]
    hrs_all = [hrs1, hrs2, hrs3]
    ilong = np.argmax([len(i) for i in hrs_all])
    # grab T data
    AERI_T = AERI_dic[d]["temperature"][:, iz]
    AERIrLID_T = AERIrLID_dic[d]["temperature"][:, iz]
    AERIvDIAL_T = AERIvDIAL_dic[d]["temperature"][:, iz]
    # grab w data
    AERI_w = AERI_dic[d]["waterVapor"][:, iz]
    AERIrLID_w = AERIrLID_dic[d]["temperature"][:, iz]
    AERIvDIAL_w = AERIvDIAL_dic[d]["temperature"][:, iz]
    # store in a dictionary for looping access
    T_dic = {0: AERI_T,
             1: AERIrLID_T,
             2: AERIvDIAL_T}
    w_dic = {0: AERI_w,
             1: AERIrLID_w,
             2: AERIvDIAL_w}
    # interpolate in time at each level
    hrs_interp = hrs_all[ilong]
    T_dic_interp = {}
    w_dic_interp = {}
    for i in range(3):
        if i == ilong:
            T_dic_interp[i] = T_dic[i]
            w_dic_interp[i] = w_dic[i]
        else:
            T_new = np.full((len(hrs_interp), len(iz)), np.nan)
            w_new = np.full((len(hrs_interp), len(iz)), np.nan)
            for j in range(len(iz)):
                fT = interp1d(hrs_all[i], T_dic[i][:, j], fill_value="extrapolate")
                fnewT = fT(hrs_interp.data)
                T_new[:, j] = fnewT
                fw = interp1d(hrs_all[i], w_dic[i][:, j], fill_value="extrapolate")
                fneww = fw(hrs_interp.data)
                w_new[:, j] = fneww
            T_dic_interp[i] = T_new
            w_dic_interp[i] = w_new

    # calculate differences
    AERI_LID_T_diff = (T_dic_interp[0] - T_dic_interp[1]).transpose()
    AERI_DIAL_T_diff = (T_dic_interp[0] - T_dic_interp[2]).transpose()
    AERI_LID_w_diff = (w_dic_interp[0] - w_dic_interp[1]).transpose()
    AERI_DIAL_w_diff = (w_dic_interp[0] - w_dic_interp[2]).transpose()
    # plot
    # Temperature
    fig2, ax2 = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(16, 8))
    # AERI - AERIrLID
    cfax21 = ax2[0].pcolormesh(hrs_interp, z[iz], AERI_LID_T_diff, 
        cmap=cmocean.cm.balance,vmin=-5,vmax=5)
    cbar21 = plt.colorbar(cfax21, ax=ax2[0])
    cbar21.ax.set_ylabel("$\Delta T$ [$^\circ$C]")
    ax2[0].set_ylabel("Altitude [km AGL]")
    ax2[0].set_ylim([0, 4])
    ax2[0].yaxis.set_major_locator(MultipleLocator(1))
    ax2[0].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax2[0].set_title(f"2017-05-{d} AERI - Raman")
    cfax21.set_edgecolor("face")

    # AERI - AERIvDIAL
    cfax22 = ax2[1].pcolormesh(hrs_interp, z[iz], AERI_DIAL_T_diff,
        cmap=cmocean.cm.balance,vmin=-5,vmax=5)
    cbar22 = plt.colorbar(cfax22, ax=ax2[1])
    cbar22.ax.set_ylabel("$\Delta T$ [$^\circ$C]")
    ax2[1].set_xlabel("Hour [UTC]")
    ax2[1].set_ylabel("Altitude [km AGL]")
    ax2[1].set_xlim([0, 24])
    ax2[1].set_ylim([0, 4])
    ax2[1].xaxis.set_major_locator(MultipleLocator(3))
    ax2[1].xaxis.set_minor_locator(MultipleLocator(1))
    ax2[1].yaxis.set_major_locator(MultipleLocator(1))
    ax2[1].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax2[1].set_title("AERI - vDIAL")
    cfax22.set_edgecolor("face")

    fig2.tight_layout()
    fsave2 = f"201705{d}_Temperature_diff.pdf"
    print(f"Saving figure: {fsave2}")
    fig2.savefig(os.path.join(figpath, fsave2), dpi=150, fmt="pdf")
    plt.close(fig2)

    # Water Vapor
    fig3, ax3 = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(16, 8))
    # AERI - AERIrLID
    cfax31 = ax3[0].pcolormesh(hrs_interp, z[iz], AERI_LID_w_diff, 
        cmap=cmocean.cm.balance,vmin=-10,vmax=10)
    cbar31 = plt.colorbar(cfax31, ax=ax3[0])
    cbar31.ax.set_ylabel("$\Delta w$ [g Kg$^{-1}$]")
    ax3[0].set_ylabel("Altitude [km AGL]")
    ax3[0].set_ylim([0, 4])
    ax3[0].yaxis.set_major_locator(MultipleLocator(1))
    ax3[0].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax3[0].set_title(f"2017-05-{d} AERI - Raman")
    cfax31.set_edgecolor("face")

    # AERI - AERIvDIAL
    cfax32 = ax3[1].pcolormesh(hrs_interp, z[iz], AERI_DIAL_w_diff,
        cmap=cmocean.cm.balance,vmin=-10,vmax=10)
    cbar32 = plt.colorbar(cfax32, ax=ax3[1])
    cbar32.ax.set_ylabel("$\Delta w$ [g Kg$^{-1}$]")
    ax3[1].set_xlabel("Hour [UTC]")
    ax3[1].set_ylabel("Altitude [km AGL]")
    ax3[1].set_xlim([0, 24])
    ax3[1].set_ylim([0, 4])
    ax3[1].xaxis.set_major_locator(MultipleLocator(3))
    ax3[1].xaxis.set_minor_locator(MultipleLocator(1))
    ax3[1].yaxis.set_major_locator(MultipleLocator(1))
    ax3[1].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax3[1].set_title("AERI - vDIAL")
    cfax32.set_edgecolor("face")

    fig3.tight_layout()
    fsave3 = f"201705{d}_waterVapor_diff.pdf"
    print(f"Saving figure: {fsave3}")
    fig3.savefig(os.path.join(figpath, fsave3), dpi=150, fmt="pdf")
    plt.close(fig3)