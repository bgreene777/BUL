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

# ---------------- #
# Data directories #
# ---------------- #
# home directory
home = os.path.expanduser("~")
# main data directory
BUL = os.path.join(home, "Documents", "Data", "BUL")

# --------------- #
# Load vDial Data #
# --------------- #
# AERIonly
f_AERI = glob(os.path.join(BUL, "vDial", "AERIonly", "*.2017052*.cdf"))
f_AERI = f_AERI[:4]
AERI_dic = {}
d = 20
# loop through these files and load all vars and ncattrs into dictionary
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
f_AERIrLID = glob(os.path.join(BUL, "vDial", "AERIrLID", "*.2017052*.cdf"))
f_AERIrLID = f_AERIrLID[:4]
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
f_AERIvDIAL = glob(os.path.join(BUL, "vDial", "AERIvDIAL", "*.2017052*.cdf"))
f_AERIvDIAL = f_AERIvDIAL[:4]
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
