'''
Load netcdf files, extract useful vars, interpolate to common time grid,
save dictionary to pickle file
Dates analyzed: 20170516 -- 20170612

Author: Brian Greene
University of Oklahoma
Last update: April 2020
'''
# Python Packages
import os
import pickle

# Installed packages
import netCDF4
import numpy as np
from datetime import datetime

from glob import glob
from scipy.interpolate import interp1d

# ---------------- #
# Data directories #
# ---------------- #
# home directory
home = os.path.expanduser("~")
# main data directory
BUL = os.path.join(home, "Documents", "Data", "BUL")
# save directory
pkl = os.path.join(os.getcwd(), "pkl")

# --------- #
# Load Data #
# --------- #
# AERIonly
f_AERI = glob(os.path.join(BUL, "WaterVaporDial", "AERIonly", "*.20170[516-609]*.cdf"))
# sort files
date_AERI = [int(ff.split(os.sep)[-1].split(".")[-3]) for ff in f_AERI]
id_AERI = np.argsort(date_AERI)
AERI_dic = {}
# loop through these files and load all vars and ncattrs into dictionary
print("---AERI Only---")
for d, f in enumerate(np.asarray(f_AERI)[id_AERI]):
    AERI_dic[d] = {}
    print(f"Reading file: {f.split(os.sep)[-1]}")
    df = netCDF4.Dataset(f, "r")
    for var in df.variables:
        AERI_dic[d][var] = df.variables[var][:]
    for attr in df.ncattrs():
        AERI_dic[d][attr] = df.getncattr(attr)
    df.close()

#AERIrLID
print("---AERI with Raman Lidar---")
f_AERIrLID = glob(os.path.join(BUL, "WaterVaporDial", "AERIrLID", "*.20170[516-609]*.cdf"))
date_AERIrLID = [int(ff.split(os.sep)[-1].split(".")[-3]) for ff in f_AERIrLID]
id_AERIrLID = np.argsort(date_AERIrLID)
AERIrLID_dic = {}
# loop through these files and load all vars and ncattrs into dictionary
for d, f in enumerate(np.asarray(f_AERIrLID)[id_AERIrLID]):
    AERIrLID_dic[d] = {}
    print(f"Reading file: {f.split(os.sep)[-1]}")
    df = netCDF4.Dataset(f, "r")
    for var in df.variables:
        AERIrLID_dic[d][var] = df.variables[var][:]
    for attr in df.ncattrs():
        AERIrLID_dic[d][attr] = df.getncattr(attr)
    df.close()

#AERIvDIAL
print("---AERI with Water Vapor Dial---")
f_AERIvDIAL = glob(os.path.join(BUL, "WaterVaporDial", "AERIvDIAL", "*.20170[516-609]*.cdf"))
date_AERIvDIAL = [int(ff.split(os.sep)[-1].split(".")[-3]) for ff in f_AERIvDIAL]
id_AERIvDIAL = np.argsort(date_AERIvDIAL)
AERIvDIAL_dic = {}
d = 0
# loop through these files and load all vars and ncattrs into dictionary
for d, f in enumerate(np.asarray(f_AERIvDIAL)[id_AERIvDIAL]):
    AERIvDIAL_dic[d] = {}
    print(f"Reading file: {f.split(os.sep)[-1]}")
    df = netCDF4.Dataset(f, "r")
    for var in df.variables:
        AERIvDIAL_dic[d][var] = df.variables[var][:]
    for attr in df.ncattrs():
        AERIvDIAL_dic[d][attr] = df.getncattr(attr)
    df.close()

# ---------------------- #
# Loop through all files #
# ---------------------- #
# want to create single 1d arrays for time, altitude, cbh plus
# 2d arrays for temperature and wvmr per retrieval
T_all = {}
w_all = {}
t_cbh_all = np.array([])
cbh_all = np.array([])
for d in range(len(f_AERI)):
    # make sure comparing same days
    dt_A = datetime.utcfromtimestamp(AERI_dic[d]["base_time"])
    dt_AR = datetime.utcfromtimestamp(AERIrLID_dic[d]["base_time"])
    dt_AD = datetime.utcfromtimestamp(AERIvDIAL_dic[d]["base_time"])
    print(f"AERI: {dt_A}, AERI+Raman: {dt_AR}, AERI+Dial: {dt_AD}")
    # grab data
    # heights below 4 km
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
    AERIrLID_w = AERIrLID_dic[d]["waterVapor"][:, iz]
    AERIvDIAL_w = AERIvDIAL_dic[d]["waterVapor"][:, iz]

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

    # if first timestep, initialize large arrays with current data
    if d == 0:
        hrs_run = np.full_like(hrs_interp, np.nan, dtype="float")
        hrs_run[:] = hrs_interp
        for i in range(3):
            T_all[i] = np.full_like(T_dic_interp[0], np.nan, dtype="float")
            T_all[i][:] = T_dic_interp[i]
            w_all[i] = np.full_like(w_dic_interp[0], np.nan, dtype="float")
            w_all[i][:] = w_dic_interp[i]
    # otherwise, concatenate to end of big arrays
    else:
        # keep time array as just days since first datapoint
        t_this = hrs_interp + (d * 24.)
        hrs_run = np.concatenate((hrs_run, t_this), axis=0)
        for i in range(3):
            T_all[i] = np.concatenate((T_all[i], T_dic_interp[i]), axis=0)
            w_all[i] = np.concatenate((w_all[i], w_dic_interp[i]), axis=0)

    # keep track of cbh
    i_cbh = np.where(AERI_dic[d]["cbh_flag"] == 2)[0]
    cbh = AERI_dic[d]["cbh"][i_cbh]
    cbh_all = np.concatenate((cbh_all, cbh), axis=0)
    t_cbh = AERI_dic[d]["hour"][i_cbh] + (d * 24.)
    t_cbh_all = np.concatenate((t_cbh_all, t_cbh), axis=0)

# -------------------------------------------- #
# Now can save pickle files for each retrieval #
# -------------------------------------------- #
# create dictionaries of data to save and write to pickle file
AERI_save = {"height": z[iz],
             "hours": hrs_run,
             "temperature": T_all[0],
             "wvmr": w_all[0],
             "cbh": cbh_all,
             "hours_cbh": t_cbh_all}
# with open(os.path.join(pkl, "aeri.pickle"), "wb") as fn:
#     pickle.dump(AERI_save, fn)

AERIrLID_save = {"height": z[iz],
                 "hours": hrs_run,
                 "temperature": T_all[1],
                 "wvmr": w_all[1],
                 "cbh": cbh_all,
                 "hours_cbh": t_cbh_all}
# with open(os.path.join(pkl, "aeri_rlid.pickle"), "wb") as fn:
#     pickle.dump(AERIrLID_save, fn)

AERIvDIAL_save = {"height": z[iz],
                  "hours": hrs_run,
                  "temperature": T_all[2],
                  "wvmr": w_all[2],
                  "cbh": cbh_all,
                  "hours_cbh": t_cbh_all}
# with open(os.path.join(pkl, "aeri_vdial.pickle"), "wb") as fn:
#     pickle.dump(AERIvDIAL_save, fn)
