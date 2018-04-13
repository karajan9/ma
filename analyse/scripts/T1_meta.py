# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.optimize import curve_fit
import glob
import os
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath("/home/karajan/uni/master/analyse/scripts"))
from nmr_lib import *


# %%

# Zürn Daten
fn_extern = "/home/karajan/uni/master/ma/analyse/data/crn/data/T1/t1-zuern.data"
zdata = np.loadtxt(fn_extern)
ztemp = zdata[:, 0]
zt1 = zdata[:, 1]


# OBI
home_dir = "/home/karajan/uni/master/ma/analyse/data/crn/data/T1/"
t1dir = [
    home_dir + "t1_230K_280K.data",
    home_dir + "t1_280K_290K.data",
    home_dir + "t1_300K_310K.data",
    home_dir + "t1_310K_355K.data",
    home_dir + "t1_342K_380K.data",
    home_dir + "t1_360K_440K_170807.data",
    home_dir + "t1_270K_330K.data",
    home_dir + "t1_305K_345K.data",
    home_dir + "t1_305K_325K.data",
]
otemps = np.empty(0)
ot1 = np.empty(0)
ot1err = np.empty(0)
obeta = np.empty(0)
obetaerr = np.empty(0)
for file_name in t1dir:
    data = np.loadtxt(file_name)

    otemps = np.hstack([otemps, data[:, 1]])
    ot1 = np.hstack([ot1, data[:, 3]])
    ot1err = np.hstack([ot1err, data[:, 4]])
    obeta = np.hstack([obeta, data[:, 5]])
    obetaerr = np.hstack([obetaerr, data[:, 6]])

cutofftemp = 390
obeta = obeta[otemps < cutofftemp]
obetaerr = obetaerr[otemps < cutofftemp]
obetatemps = otemps[otemps < cutofftemp]


# Bruker
bruker_name = "/home/karajan/uni/master/ma/analyse/data/crn/data/T1/bruker_t1.data"
bruker_data = np.loadtxt(bruker_name)
bruker_temp = bruker_data[:, 1]
bruker_t1 = bruker_data[:, 3]
bruker_t1_err = bruker_data[:, 4]
bruker_beta = bruker_data[:, 5]
bruker_beta_err = bruker_data[:, 6]


# %%
# gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
# plt.figure(figsize=(6, 6))
f, (pltbeta, pltt1) = plt.subplots(2, 1, sharex=True,
                                   gridspec_kw={"height_ratios": [1, 4]})
f.subplots_adjust(hspace=0.1)
f.set_size_inches(6, 5)

# pltbeta = plt.subplot(gs[0])
# pltt1 = plt.subplot(gs[1])


pltt1.scatter(
    ztemp,
    1 / zt1,
    color="tab:orange",
    # facecolors="none",
    marker="x",
    label="Zürn")

pltt1.errorbar(
    otemps,
    1 / ot1,
    yerr=ot1err / ot1**2,
    linestyle="None",
    color="tab:blue",
    marker="o",
    label="OBI")

pltt1.errorbar(
    bruker_temp,
    1 / bruker_t1,
    yerr=bruker_t1_err / bruker_t1**2,
    linestyle="None",
    color="tab:red",
    marker="D",
    label="Bruker")

pltt1.set_yscale("log")
pltt1.set_xlim(220, 445)
pltt1.set_ylim(30, 10**6)
pltt1.tick_params(which="both", direction="in", top=True, right=True)
pltt1.set_xlabel("Temperatur [K]")
pltt1.set_ylabel("1/$T_1$ [1/s]")
pltt1.legend(loc=2)


pltbeta.errorbar(
    obetatemps,
    obeta,
    yerr=obetaerr,
    linestyle="None",
    color="tab:blue",
    markersize=4,
    marker="o")

pltbeta.errorbar(
    bruker_temp,
    bruker_beta,
    yerr=bruker_beta_err,
    linestyle="None",
    color="tab:red",
    markersize=4,
    marker="D")

pltbeta.set_xlim(220, 445)
pltbeta.set_ylim(0.4, 1.1)
pltbeta.tick_params(which="both", direction="in", top=True, right=True,
                    labelbottom=False, labeltop=True)
pltbeta.set_ylabel("$\\beta$")


save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/t1")
