# %%
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import gridspec
import matplotlib
from scipy.optimize import curve_fit
import glob
import os
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(
    os.path.abspath("/home/karajan/uni/master/ma/analyse/scripts"))
from nmr_lib import *


# %%

# Zürn Daten
home_dir = "/home/karajan/uni/master/ma/analyse"
fn_zuernt1 = home_dir + "/data/crn/data/T1/t1-zuern.data"
zdatat1 = np.loadtxt(fn_zuernt1)
ztempt1 = zdatat1[:, 0]
zt1 = zdatat1[:, 1]

fn_zuernbeta = home_dir + "/data/crn/data/T1/t1beta-zuern.data"
zdatabeta = np.loadtxt(fn_zuernbeta)
ztempbeta = zdatabeta[:, 0]
zbeta = zdatabeta[:, 1]


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


fn_t1betafest = home_dir + "T1_170807_betafest.data"
data_t1_bf = np.loadtxt(fn_t1betafest)

otemps_bf = data_t1_bf[:, 1]
ot1_bf = data_t1_bf[:, 3]
ot1err_bf = data_t1_bf[:, 4]
obeta_bf = data_t1_bf[:, 5]
obetaerr_bf = data_t1_bf[:, 6]


cutofftemp = 390
ot1 = ot1[otemps < cutofftemp]
ot1err = ot1err[otemps < cutofftemp]
obeta = obeta[otemps < cutofftemp]
obetaerr = obetaerr[otemps < cutofftemp]
otemps = otemps[otemps < cutofftemp]

ot1_bf = ot1_bf[otemps_bf > cutofftemp]
ot1err_bf = ot1err_bf[otemps_bf > cutofftemp]
obeta_bf = obeta_bf[otemps_bf > cutofftemp]
obetaerr_bf = obetaerr_bf[otemps_bf > cutofftemp]
otemps_bf = otemps_bf[otemps_bf > cutofftemp]


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


# plt.rc('text', usetex=True)
pltt1.errorbar(
    otemps,
    1 / ot1,
    yerr=ot1err / ot1**2,
    linestyle="None",
    color="tab:blue",
    marker="o",
    label="$f_{L}(\\mathrm{OBI}) = 97.2$ MHz")

pltt1.errorbar(
    otemps_bf,
    1 / ot1_bf,
    yerr=ot1err_bf / ot1_bf**2,
    linestyle="None",
    color="tab:blue",
    # label="OBI",
    marker="*",
)

pltt1.errorbar(
    bruker_temp,
    1 / bruker_t1,
    yerr=bruker_t1_err / bruker_t1**2,
    linestyle="None",
    color="tab:red",
    marker="D",
    label="$f_{L}(\\mathrm{Bruker}) = 131.0$ MHz")

pltt1.errorbar(
    ztempt1,
    1 / zt1,
    color="tab:orange",
    linestyle="None",
    # markersize=7,
    # facecolors="none",
    marker="^",
    label="$f_{L}(\\mathrm{Zürn}) = 85.7$ MHz")


pltt1.set_yscale("log")
pltt1.set_xlim(220, 445)
pltt1.set_ylim(30, 10**6)
pltt1.tick_params(which="both", direction="in", top=True, right=True)
pltt1.set_xlabel("Temperatur [K]")
pltt1.set_ylabel("1/$T_1$ [1/s]")
pltt1.legend(loc=2)




pltbeta.errorbar(
    otemps,
    obeta,
    yerr=obetaerr,
    linestyle="None",
    color="tab:blue",
    markersize=4,
    marker="o")

pltbeta.errorbar(
    otemps_bf,
    obeta_bf,
    yerr=obetaerr_bf,
    linestyle="None",
    color="tab:blue",
    markersize=4,
    marker="*")

pltbeta.errorbar(
    bruker_temp,
    bruker_beta,
    yerr=bruker_beta_err,
    linestyle="None",
    color="tab:red",
    markersize=4,
    marker="D")

pltbeta.errorbar(
    ztempbeta,
    zbeta,
    linestyle="None",
    color="tab:orange",
    markersize=4,
    marker="^")


pltbeta.set_xlim(220, 445)
pltbeta.set_ylim(0.4, 1.1)
pltbeta.tick_params(which="both", direction="in", top=True, right=True,
                    labelbottom=False, labeltop=True)
pltbeta.set_ylabel("$\\beta$")


save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/t1")
