'''
Produce composite difference versus height plots comparing AERIoe retrievals
Dates analyzed: 20170516 -- 20170612

Author: Brian Greene
University of Oklahoma
Last update: April 2020
'''
# Python Packages
import os
from datetime import datetime, timedelta
import warnings
import pickle
import gzip

# Installed packages
import netCDF4
import numpy as np
import cmocean

from glob import glob
from datetime import datetime
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
# piclke file data (same directory as this file)
pkl = os.path.join(os.getcwd(), "pkl")
# figure save directory
figpath = os.path.join(BUL, "Figures", "composite")
if not os.path.exists(figpath):
    os.mkdir(figpath)
# close all figures
plt.close("all")

# -------------------------------------------------------------- #
# Load 3 pickle files as dictionaries sorted by date of the data #
# -------------------------------------------------------------- #
# define reading function
def read_pickle_gz(filepath):
    print(f"Reading file: {filepath.split(os.sep)[-1]}")
    with gzip.open(filepath, "rb") as f:
        data = pickle.load(f)
    return data

def read_pickle(filepath):
    print(f"Reading file: {filepath.split(os.sep)[-1]}")
    with open(filepath, "rb") as f:
        data = pickle.load(f)
    return data

# load data
f_AERI = os.path.join(pkl, "aeri.pickle")
data_a = read_pickle(f_AERI)
f_AERIrLID = os.path.join(pkl, "aeri_rlid.pickle")
data_ar = read_pickle(f_AERIrLID)
f_AERIvDIAL = os.path.join(pkl, "aeri_vdial.pickle")
data_av = read_pickle(f_AERIvDIAL)

# -------------------- #
# Calculate statistics #
# -------------------- #
# temperature
# raman
T_diff_ar = data_a["temperature"] - data_ar["temperature"]
T_med_ar = np.nanmedian(T_diff_ar, axis=0)
T_q1_ar = np.nanpercentile(T_diff_ar, 25., axis=0)
T_q3_ar = np.nanpercentile(T_diff_ar, 75., axis=0)
# wv dial
T_diff_av = data_a["temperature"] - data_av["temperature"]
T_med_av = np.nanmedian(T_diff_av, axis=0)
T_q1_av = np.nanpercentile(T_diff_av, 25., axis=0)
T_q3_av = np.nanpercentile(T_diff_av, 75., axis=0)

# wvmr
# raman
w_diff_ar = data_a["wvmr"] - data_ar["wvmr"]
w_med_ar = np.nanmedian(w_diff_ar, axis=0)
w_q1_ar = np.nanpercentile(w_diff_ar, 25., axis=0)
w_q3_ar = np.nanpercentile(w_diff_ar, 75., axis=0)
# wv dial
w_diff_av = data_a["wvmr"] - data_av["wvmr"]
w_med_av = np.nanmedian(w_diff_av, axis=0)
w_q1_av = np.nanpercentile(w_diff_av, 25., axis=0)
w_q3_av = np.nanpercentile(w_diff_av, 75., axis=0)

# ---- #
# Plot #
# ---- #
rc('font',weight='normal',size=20,family='serif',serif='Computer Modern Roman')
rc('text',usetex='True')
z = data_a["height"]

fig1, ax1 = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 8))
# temperature differences
# raman
ax1[0].plot(T_med_ar, z, "-b", linewidth=3.)
ax1[0].fill_betweenx(z, T_q1_ar, T_q3_ar, alpha=0.3, color="b")
# wv dial
ax1[0].plot(T_med_av, z, "-r", linewidth=3.)
ax1[0].fill_betweenx(z, T_q1_av, T_q3_av, alpha=0.3, color="r")
# subplot setup
ax1[0].axvline(0., linewidth=2., color="k", linestyle="--")
ax1[0].grid()
ax1[0].set_xlabel(r"$T_{AERI} - (T_{AERIrLID}, T_{AERIvDIAL})$ [$^\circ$C]", fontsize=20)
ax1[0].set_ylabel("Height [m AGL]", fontsize=20)
ax1[0].set_xlim([-0.4, 0.6])
ax1[0].set_ylim([0., 3.7])

# wvmr differences
ax1[1].plot(w_med_ar, z, "-b", linewidth=3)
ax1[1].fill_betweenx(z, w_q1_ar, w_q3_ar, alpha=0.3, color="b")
# wv dial
ax1[1].plot(w_med_av, z, "-r", linewidth=3)
ax1[1].fill_betweenx(z, w_q1_av, w_q3_av, alpha=0.3, color="r")
# subplot setup
ax1[1].axvline(0., linewidth=2., color="k", linestyle="--")
ax1[1].grid()
ax1[1].set_xlabel(r"$w_{AERI} - (w_{AERIrLID}, w_{AERIvDIAL})$ [g kg$^{-1}$]", fontsize=20)
ax1[1].set_xlim([-1.5, 1.5])

#
# Plot all diffs as xs
#
# vdial
fig2, ax2 = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 8))
# temperature differences
for i in range(len(z)):
    ax2[0].plot(T_diff_av[:, i], z[i]*np.ones(len(T_diff_av[:, i])), "kx", alpha=0.3)
ax2[0].plot(T_med_ar, z, "ro")
ax2[0].plot(np.nanmean(T_diff_av, axis=0), z, "bx")
ax2[0].grid()
ax2[0].set_ylim([0., 1.])
ax2[0].set_xlim([-6., 6.])
ax2[0].set_ylabel("Height [km AGL]")
ax2[0].set_xlabel("AERIonly T - AERIvDIAL T [$^\circ$C]")

# wvmr diff
for i in range(len(z)):
    ax2[1].plot(w_diff_av[:, i], z[i]*np.ones(len(w_diff_av[:, i])), "kx", alpha=0.3)
ax2[1].plot(w_med_ar, z, "ro")
ax2[1].plot(np.nanmean(w_diff_av, axis=0), z, "bx")
ax2[1].grid()
ax2[1].set_xlim([-6., 6.])
ax2[1].set_xlabel("AERIonly w - AERIvDIAL w [g kg$^{-1}]")

plt.show()