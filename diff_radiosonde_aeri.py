'''
Produce composite difference versus height plots comparing AERIoe retrievals
to radiosondes launched at ARM SGP
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
figpath = os.path.join(BUL, "Figures", "radiosondes")
if not os.path.exists(figpath):
    os.mkdir(figpath)
# close all figures
plt.close("all")

# ---------------------------------------------------------- #
# Load 3 pickle files as dictionaries with time/height grids #
# ---------------------------------------------------------- #
# define reading function
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

# ----------------------------------- #
# Loop through all sondes and compare #
# ----------------------------------- #
f_sond = glob(os.path.join(BUL, "Radiosonde", "*.cdf"))
# sort by date and time on filename
dt_sond = [datetime.strptime("".join(ff.split(os.sep)[-1].split(".")[-3:-1]), 
    "%Y%m%d%H%M%S") for ff in f_sond]
i_sond = np.argsort(dt_sond)

# initialize lists for T, Td, alt, time
T_all = []
w_all = []
alt_all = []
t_all = []
# define some datetime object constants
dt_0 = datetime(2017, 5, 16)
delta = timedelta(seconds=1)

# define function to convert Td and p into wvmr
def dewpoint_to_wvmr(Td, p):
    '''
    input: Td in C, p in hPa
    output: wvmr in g/kg
    use following eqns
    e = 6.112 * exp(17.67*Td / (Td + 243.5))
    w = 0.622 * e / (p - e)
    '''
    e = 6.112 * np.exp(17.67 * Td / (Td + 243.5))
    w = 0.622 * e / (p - e)
    return 1000. * w

# define function to locate nearest aeri point
def find_nearest_timestep(t_comp, t_ref):
    '''
    input: t_comp = single timestep to compare against all reference times
    input: t_ref = array of timesteps to search through
    output: index of reference timestep in t_ref closest to t_comp
    output: how close (in same units of time) to this point
    '''
    return np.argmin(np.abs([(t_comp - it) for it in t_ref]))

# begin loop to load data
for i, f in enumerate(np.array(f_sond)[i_sond]):
    # read
    print(f"Reading file: {f.split(os.sep)[-1]}")
    df = netCDF4.Dataset(f, "r")
    # only grab data from below 4 km to compare
    i_4km = np.where(df.variables["alt"][:].data < 4000.)[0]
    # append to lists
    T_all.append(df.variables["tdry"][i_4km].data)
    alt_all.append(df.variables["alt"][i_4km].data)
    # calculate wvmr from Td and p
    w = dewpoint_to_wvmr(df.variables["dp"][i_4km].data, 
                         df.variables["pres"][i_4km].data)
    w_all.append(w)
    # convert time data to be consistent notation as AERI
    # each file records time as seconds since 0Z on same day
    # create datetime objects for each timestep
    dt_this = np.array(dt_sond)[i_sond][i]
    dt_base = datetime(dt_this.year, dt_this.month, dt_this.day)
    t = df.variables["time"][i_4km].data
    dt_new = np.array([(dt_base + (n * delta)) for n in t])
    # now convert to hours since first day at 0Z to be same as AERI
    t_new = np.array([((it - dt_0).total_seconds()/3600.) for it in dt_new])
    # all done; can append to t_all
    t_all.append(t_new)
    # close file
    df.close()

# define new altitude grid in m
# convert from AGL to MSL for comparison with sondes
z_a = (data_a["height"].data * 1000.) + 237.43
nz = len(z_a)
# initialize arrays to fill with gridded interpolated data
n_sond = len(f_sond)
T_grid = np.full((n_sond, nz), np.nan, dtype=float)
w_grid = np.full((n_sond, nz), np.nan, dtype=float)
t_grid = np.full((n_sond, nz), np.nan, dtype=float)
# initialize arrays to store closest AERI data
# aeri
T_close_a = np.full((n_sond, nz), np.nan, dtype=float)
w_close_a = np.full((n_sond, nz), np.nan, dtype=float)
# aeri + raman
T_close_ar = np.full((n_sond, nz), np.nan, dtype=float)
w_close_ar = np.full((n_sond, nz), np.nan, dtype=float)
# aeri + wv dial
T_close_av = np.full((n_sond, nz), np.nan, dtype=float)
w_close_av = np.full((n_sond, nz), np.nan, dtype=float)
# begin loop to calculate comparisons
for i in range(n_sond):
    # want to get radiosonde data in same gridded format as AERI
    # interpolate to z_a
    # temperature
    fT = interp1d(alt_all[i], T_all[i], fill_value=np.nan, bounds_error=False)
    fnewT = fT(z_a)
    T_grid[i, :] = fnewT
    # wvmr
    fw = interp1d(alt_all[i], w_all[i], fill_value=np.nan, bounds_error=False)
    fneww = fw(z_a)
    w_grid[i, :] = fneww
    # time
    ft = interp1d(alt_all[i], t_all[i], fill_value=np.nan, bounds_error=False)
    fnewt = ft(z_a)
    t_grid[i, :] = fnewt
    # now can loop through each gridded timestep to compare with AERIoe
    for j, jt in enumerate(fnewt):
        i_close = find_nearest_timestep(jt, data_a["hours"])
        T_close_a[i, j] = data_a["temperature"][i_close, j]
        w_close_a[i, j] = data_a["wvmr"][i_close, j]
        T_close_ar[i, j] = data_ar["temperature"][i_close, j]
        w_close_ar[i, j] = data_ar["wvmr"][i_close, j]
        T_close_av[i, j] = data_av["temperature"][i_close, j]
        w_close_av[i, j] = data_av["wvmr"][i_close, j]

    # # better way: find closest timestep of lowest available data point
    # # and find out how close to the next timestep in indices, then can
    # # more efficiently choose timestamps to compare
    # # first point always at iz = 7
    # i_close, t_close = find_nearest_timestep(fnewt[7], data_a["hours"])
    # # find time halfway between i_close and i_close+1
    # t_mid = np.mean(data_a["hours"][i_close:i_close+2])

# calculate differences and statistics
def RMSD(obs, ref):
    return np.sqrt(np.nanmean((obs - ref)**2.))
# AERIonly
# temperature
T_rmsd_a = RMSD(T_grid, T_close_a)
T_diff_a = T_grid - T_close_a
T_med_a = np.nanmedian(T_diff_a, axis=0)
T_q1_a = np.nanpercentile(T_diff_a, 25., axis=0)
T_q3_a = np.nanpercentile(T_diff_a, 75., axis=0)
# wvmr
w_rmsd_a = RMSD(w_grid, w_close_a)
w_diff_a = w_grid - w_close_a
w_med_a = np.nanmedian(w_diff_a, axis=0)
w_q1_a = np.nanpercentile(w_diff_a, 25., axis=0)
w_q3_a = np.nanpercentile(w_diff_a, 75., axis=0)
# AERI + raman
# temperature
T_rmsd_ar = RMSD(T_grid, T_close_ar)
T_diff_ar = T_grid - T_close_ar
T_med_ar = np.nanmedian(T_diff_ar, axis=0)
T_q1_ar = np.nanpercentile(T_diff_ar, 25., axis=0)
T_q3_ar = np.nanpercentile(T_diff_ar, 75., axis=0)
# wvmr
w_rmsd_ar = RMSD(w_grid, w_close_ar)
w_diff_ar = w_grid - w_close_ar
w_med_ar = np.nanmedian(w_diff_ar, axis=0)
w_q1_ar = np.nanpercentile(w_diff_ar, 25., axis=0)
w_q3_ar = np.nanpercentile(w_diff_ar, 75., axis=0)
# AERI + wv dial
# temperature
T_rmsd_av = RMSD(T_grid, T_close_av)
T_diff_av = T_grid - T_close_av
T_med_av = np.nanmedian(T_diff_av, axis=0)
T_q1_av = np.nanpercentile(T_diff_av, 25., axis=0)
T_q3_av = np.nanpercentile(T_diff_av, 75., axis=0)
# wvmr
w_rmsd_av = RMSD(w_grid, w_close_av)
w_diff_av = w_grid - w_close_av
w_med_av = np.nanmedian(w_diff_av, axis=0)
w_q1_av = np.nanpercentile(w_diff_av, 25., axis=0)
w_q3_av = np.nanpercentile(w_diff_av, 75., axis=0)

# calculate 2d histogram bins and edges
# temperature
T_bins = (np.arange(-10., 35.5, 0.5), np.arange(-10., 35.5, 0.5))
H_a, xe_a, ye_a = np.histogram2d(np.ravel(T_close_a),
                                 np.ravel(T_grid),
                                 bins=T_bins,
                                 density=True)
xaT, yaT = np.meshgrid(xe_a, ye_a)
H_ar, xe_ar, ye_ar = np.histogram2d(np.ravel(T_close_ar),
                                    np.ravel(T_grid),
                                    bins=T_bins,
                                    density=True)
xarT, yarT = np.meshgrid(xe_ar, ye_ar)
H_av, xe_av, ye_av = np.histogram2d(np.ravel(T_close_av),
                                    np.ravel(T_grid),
                                    bins=T_bins,
                                    density=True)
xavT, yavT = np.meshgrid(xe_av, ye_av)
# wvmr
w_bins = (np.arange(0., 18.5, 0.5), np.arange(0., 18.5, 0.5))
W_a, xew_a, yew_a = np.histogram2d(np.ravel(w_close_a),
                                   np.ravel(w_grid),
                                   bins=w_bins,
                                   density=True)
xaw, yaw = np.meshgrid(xew_a, yew_a)
W_ar, xew_ar, yew_ar = np.histogram2d(np.ravel(w_close_ar),
                                      np.ravel(w_grid),
                                      bins=w_bins,
                                      density=True)
xarw, yarw = np.meshgrid(xew_ar, yew_ar)
W_av, xew_av, yew_av = np.histogram2d(np.ravel(w_close_av),
                                      np.ravel(w_grid),
                                      bins=w_bins,
                                      density=True)
xavw, yavw = np.meshgrid(xew_av, yew_av)

# -------------------- #
# Plot AERIoe - Sondes #
# -------------------- #
print("Begin plotting...")
rc('font',weight='normal',size=20,family='serif',serif='Computer Modern Roman')
rc('text',usetex='True')
z_a_agl = z_a - 237.43
colors = [(0., 0., 0.), (0./255, 114./255., 178./255), (213./255, 94./255, 0.)]

fig1, ax1 = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 8))
# temperature
# AERI only
ax1[0].plot(T_med_a, z_a_agl, color=colors[0], linestyle="-", linewidth=3., 
    label="AERIonly")
# ax1[0].fill_betweenx(z_a_agl, T_q1_a, T_q3_a, alpha=0.3, color=colors[0])
ax1[0].plot(T_q1_a, z_a_agl, T_q3_a, z_a_agl, color=colors[0], linestyle=":", linewidth=2)
# AERI + Raman
ax1[0].plot(T_med_ar, z_a_agl, color=colors[1], linestyle="-", linewidth=3., 
    label="AERIrLID")
# ax1[0].fill_betweenx(z_a_agl, T_q1_ar, T_q3_ar, alpha=0.3, color=colors[1])
ax1[0].plot(T_q1_ar, z_a_agl, T_q3_ar, z_a_agl, color=colors[1], linestyle=":", linewidth=2)
# AERI + wv dial
ax1[0].plot(T_med_av, z_a_agl, color=colors[2], linestyle="-", linewidth=3., 
    label="AERIvDial")
# ax1[0].fill_betweenx(z_a_agl, T_q1_av, T_q3_av, alpha=0.3, color=colors[2])
ax1[0].plot(T_q1_av, z_a_agl, T_q3_av, z_a_agl, color=colors[2], linestyle=":", linewidth=2)
# setup
ax1[0].axvline(0., linewidth=2., color="k", linestyle="--")
ax1[0].grid()
ax1[0].legend(loc=(0.01, 0.75), fontsize=16)
ax1[0].set_xlabel("$T_{sonde} - T_{AERIoe}$ [$^\circ$C]")
ax1[0].set_ylabel("Altitude [m AGL]")
ax1[0].set_xlim([-2., 2.])
ax1[0].set_ylim([0., 4000.])
ax1[0].xaxis.set_minor_locator(MultipleLocator(0.25))
ax1[0].yaxis.set_minor_locator(MultipleLocator(250))
props=dict(boxstyle='square',facecolor='white',alpha=0.5)
ax1[0].text(0.03,0.95,r'\textbf{(a)}',fontsize=20,bbox=props, transform=ax1[0].transAxes)

# wvmr
# AERI only
ax1[1].plot(w_med_a, z_a_agl, color=colors[0], linestyle="-", linewidth=3.)
# ax1[1].fill_betweenx(z_a_agl, w_q1_a, w_q3_a, alpha=0.3, color=colors[0])
ax1[1].plot(w_q1_a, z_a_agl, w_q3_a, z_a_agl, color=colors[0], linestyle=":", linewidth=2)
# AERI + Raman
ax1[1].plot(w_med_ar, z_a_agl, color=colors[1], linestyle="-", linewidth=3.)
# ax1[1].fill_betweenx(z_a_agl, w_q1_ar, w_q3_ar, alpha=0.3, color=colors[1])
ax1[1].plot(w_q1_ar, z_a_agl, w_q3_ar, z_a_agl, color=colors[1], linestyle=":", linewidth=2)
# AERI + wv dial
ax1[1].plot(w_med_av, z_a_agl, color=colors[2], linestyle="-", linewidth=3.)
# ax1[1].fill_betweenx(z_a_agl, w_q1_av, w_q3_av, alpha=0.3, color=colors[2])
ax1[1].plot(w_q1_av, z_a_agl, w_q3_av, z_a_agl, color=colors[2], linestyle=":", linewidth=2)
# setup
ax1[1].axvline(0., linewidth=2., color="k", linestyle="--")
ax1[1].grid()
ax1[1].set_xlabel("$WVMR_{sonde} - WVMR_{AERIoe}$ [g kg$^{-1}$]")
ax1[1].set_xlim([-2., 2.])
ax1[1].xaxis.set_minor_locator(MultipleLocator(0.25))
ax1[1].yaxis.set_minor_locator(MultipleLocator(250))
ax1[1].text(0.03,0.95,r'\textbf{(b)}',fontsize=20,bbox=props, transform=ax1[1].transAxes)


fig1.tight_layout()
fsave1 = "diff_vs_alt_T_wvmr"
fig1.savefig(f"{os.path.join(figpath, fsave1)}.pdf", dpi=300, fmt="pdf")
fig1.savefig(f"{os.path.join(figpath, fsave1)}.png", dpi=300, fmt="png")
plt.close(fig1)

# ------------------ #
# Plot 2d histograms #
# ------------------ #
inset = False
fig2 = plt.figure(constrained_layout=True, figsize=(16, 10))
spec2 = fig2.add_gridspec(nrows=2, ncols=3)
# Temperature
# AERI only
ax21 = fig2.add_subplot(spec2[0, 0])
cfax21 = ax21.pcolormesh(xaT, yaT, H_a.transpose(), cmap=cmocean.cm.amp)
cfax21.set_edgecolor("face")
ax21.plot(T_bins[0], T_bins[1], color="k", linewidth=1, linestyle=":")
ax21.xaxis.set_minor_locator(MultipleLocator(5))
ax21.yaxis.set_minor_locator(MultipleLocator(5))
ax21.grid(which="both")
ax21.text(0.05,0.9,r'\textbf{(a)}',fontsize=20,bbox=props, transform=ax21.transAxes)
ax21.set_xlabel("$T_{AERIonly}$ [$^\circ$C]")
ax21.set_ylabel("$T_{sonde}$ [$^\circ$C]")
ax21.text(0.05, 0.78, f"RMSD: {T_rmsd_a:3.2f}$^\circ$C", fontsize=20,bbox=props,transform=ax21.transAxes)
# inset distribution
if inset:
    ax21_1 = inset_axes(ax21, width="35%", height="35%", loc=4, borderpad=1)
    ax21_1.hist(np.ravel(T_diff_a), bins=np.arange(-10., 10., 0.5), density=True)
    ax21_1.tick_params(axis="both", labelsize=10)
    ax21_1.set_title("$T_{sonde}-T_{AERI}$", fontsize=12)

# AERI + raman
ax22 = fig2.add_subplot(spec2[0, 1], sharey=ax21)
cfax22 = ax22.pcolormesh(xarT, yarT, H_ar.transpose(), cmap=cmocean.cm.amp)
cfax22.set_edgecolor("face")
ax22.tick_params(axis="y", labelleft=False)
ax22.plot(T_bins[0], T_bins[1], color="k", linewidth=1, linestyle=":")
ax22.xaxis.set_minor_locator(MultipleLocator(5))
ax22.yaxis.set_minor_locator(MultipleLocator(5))
ax22.grid(which="both")
ax22.text(0.05,0.9,r'\textbf{(b)}',fontsize=20,bbox=props, transform=ax22.transAxes)
ax22.set_xlabel("$T_{AERIrLID}$ [$^\circ$C]")
ax22.text(0.05, 0.78, f"RMSD: {T_rmsd_ar:3.2f}$^\circ$C", fontsize=20,bbox=props,transform=ax22.transAxes)
# inset distribution
if inset:
    ax22_1 = inset_axes(ax22, width="35%", height="35%", loc=4, borderpad=1)
    ax22_1.hist(np.ravel(T_diff_ar), bins=np.arange(-10., 10., 0.5), density=True)
    ax22_1.tick_params(axis="both", labelsize=10)
    ax22_1.set_title("$T_{sonde}-T_{AERIrLID}$", fontsize=12)

# AERI + wv dial
ax23 = fig2.add_subplot(spec2[0, 2], sharey=ax21)
cfax23 = ax23.pcolormesh(xavT, yavT, H_av.transpose(), cmap=cmocean.cm.amp)
cfax23.set_edgecolor("face")
ax23.tick_params(axis="y", labelleft=False)
ax23.plot(T_bins[0], T_bins[1], color="k", linewidth=1, linestyle=":")
ax23.xaxis.set_minor_locator(MultipleLocator(5))
ax23.yaxis.set_minor_locator(MultipleLocator(5))
ax23.grid(which="both")
ax23.text(0.05,0.9,r'\textbf{(c)}',fontsize=20,bbox=props, transform=ax23.transAxes)
ax23.set_xlabel("$T_{AERIvDIAL}$ [$^\circ$C]")
ax23.text(0.05, 0.78, f"RMSD: {T_rmsd_av:3.2f}$^\circ$C", fontsize=20,bbox=props,transform=ax23.transAxes)
# inset distribution
if inset:
    ax23_1 = inset_axes(ax23, width="35%", height="35%", loc=4, borderpad=1)
    ax23_1.hist(np.ravel(T_diff_av), bins=np.arange(-10., 10., 0.5), density=True)
    ax23_1.tick_params(axis="both", labelsize=10)
    ax23_1.set_title("$T_{sonde}-T_{AERIvDIAL}$", fontsize=12)

# wvmr
# AERI only
ax24 = fig2.add_subplot(spec2[1, 0])
cfax24 = ax24.pcolormesh(xaw, yaw, W_a.transpose(), cmap=cmocean.cm.rain)
cfax24.set_edgecolor("face")
ax24.plot(w_bins[0], w_bins[1], color="k", linewidth=1, linestyle=":")
ax24.xaxis.set_minor_locator(MultipleLocator(1))
ax24.yaxis.set_minor_locator(MultipleLocator(1))
ax24.grid(which="both")
ax24.text(0.05,0.9,r'\textbf{(d)}',fontsize=20,bbox=props, transform=ax24.transAxes)
ax24.set_xlabel("$WVMR_{AERIonly}$ [g kg$^{-1}$]")
ax24.set_ylabel("$WVMR_{sonde}$ [g kg$^{-1}$]")
ax24.text(0.05, 0.78, f"RMSD: {w_rmsd_a:3.2f} "+ r"g kg$^{-1}$",fontsize=20,bbox=props,transform=ax24.transAxes)
# inset distribution
if inset:
    ax24_1 = inset_axes(ax24, width="35%", height="35%", loc=4, borderpad=1)
    ax24_1.hist(np.ravel(w_diff_a), bins=np.arange(-10., 10., 0.5), density=True)
    ax24_1.tick_params(axis="both", labelsize=10)
    ax24_1.set_title("$WVMR_{sonde}-WVMR_{AERI}$", fontsize=12)

# AERI + Raman
ax25 = fig2.add_subplot(spec2[1, 1], sharey=ax24)
cfax25 = ax25.pcolormesh(xarw, yarw, W_ar.transpose(), cmap=cmocean.cm.rain)
cfax25.set_edgecolor("face")
ax25.plot(w_bins[0], w_bins[1], color="k", linewidth=1, linestyle=":")
ax25.xaxis.set_minor_locator(MultipleLocator(1))
ax25.yaxis.set_minor_locator(MultipleLocator(1))
ax25.grid(which="both")
ax25.text(0.05,0.9,r'\textbf{(e)}',fontsize=20,bbox=props, transform=ax25.transAxes)
ax25.set_xlabel("$WVMR_{AERIrLID}$ [g kg$^{-1}$]")
ax25.text(0.05, 0.78, f"RMSD: {w_rmsd_ar:3.2f} "+ r"g kg$^{-1}$",fontsize=20,bbox=props,transform=ax25.transAxes)
# inset distribution
if inset:
    ax25_1 = inset_axes(ax25, width="35%", height="35%", loc=4, borderpad=1)
    ax25_1.hist(np.ravel(w_diff_ar), bins=np.arange(-10., 10., 0.5), density=True)
    ax25_1.tick_params(axis="both", labelsize=10)
    ax25_1.set_title("$WVMR_{sonde}-WVMR_{AERIrLID}$", fontsize=12)

# AERI + wv dial
ax26 = fig2.add_subplot(spec2[1, 2], sharey=ax24)
cfax26 = ax26.pcolormesh(xavw, yavw, W_av.transpose(), cmap=cmocean.cm.rain)
cfax26.set_edgecolor("face")
ax26.plot(w_bins[0], w_bins[1], color="k", linewidth=1, linestyle=":")
ax26.xaxis.set_minor_locator(MultipleLocator(1))
ax26.yaxis.set_minor_locator(MultipleLocator(1))
ax26.grid(which="both")
ax26.text(0.05,0.9,r'\textbf{(f)}',fontsize=20,bbox=props, transform=ax26.transAxes)
ax26.set_xlabel("$WVMR_{AERIvDIAL}$ [g kg$^{-1}$]")
ax26.text(0.05, 0.78, f"RMSD: {w_rmsd_av:3.2f} "+ r"g kg$^{-1}$",fontsize=20,bbox=props,transform=ax26.transAxes)
# inset distribution
if inset:
    ax26_1 = inset_axes(ax26, width="35%", height="35%", loc=4, borderpad=1)
    ax26_1.hist(np.ravel(w_diff_av), bins=np.arange(-10., 10., 0.5), density=True)
    ax26_1.tick_params(axis="both", labelsize=10)
    ax26_1.set_title("$WVMR_{sonde}-WVMR_{AERIvDIAL}$", fontsize=12)

# save
fsave2 = "2d_hist_T_wvmr_all"
fig2.savefig(f"{os.path.join(figpath, fsave2)}.pdf", fmt="pdf")
fig2.savefig(f"{os.path.join(figpath, fsave2)}.png", dpi=300, fmt="png")
plt.close(fig2)

plt.close("all")