# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath("/home/karajan/uni/master/analyse/scripts"))
from nmr_lib import *


# %%

# OBI
home_dir = "/home/karajan/uni/master/ma/analyse/data/crn/data/T2/"
fn = [
    home_dir + "t2_230K_280K.data",
    home_dir + "t2_300K_310K.data",
    home_dir + "t2_360K_440K.data",
    home_dir + "t2_342K_380K.data",
    home_dir + "t2_270K_330K.data",
    home_dir + "t2_305K_345K.data",
    home_dir + "t2_305K_325K.data",
]
otemps = np.empty(0)
ot2 = np.empty(0)
ot2err = np.empty(0)
obeta = np.empty(0)
obetaerr = np.empty(0)
for file_name in fn:
    data = np.loadtxt(file_name)

    otemps = np.hstack([otemps, data[:, 1]])
    ot2 = np.hstack([ot2, data[:, 3]])
    ot2err = np.hstack([ot2err, data[:, 4]])
    obeta = np.hstack([obeta, data[:, 5]])
    obetaerr = np.hstack([obetaerr, data[:, 6]])


fn_t2betafest = home_dir + "T2_170807_betafest.data"
data_t2_bf = np.loadtxt(fn_t2betafest)

otemps_bf = data_t2_bf[:, 1]
ot2_bf = data_t2_bf[:, 3]
ot2err_bf = data_t2_bf[:, 4]
obeta_bf = data_t2_bf[:, 5]
obetaerr_bf = data_t2_bf[:, 6]


cutofftemp = 390
# ot2 = ot2[otemps < cutofftemp]
# ot2err = ot2err[otemps < cutofftemp]
# obeta = obeta[otemps < cutofftemp]
# obetaerr = obetaerr[otemps < cutofftemp]
# otemps = otemps[otemps < cutofftemp]

ot2_bf = ot2_bf[otemps_bf > cutofftemp]
ot2err_bf = ot2err_bf[otemps_bf > cutofftemp]
obeta_bf = obeta_bf[otemps_bf > cutofftemp]
obetaerr_bf = obetaerr_bf[otemps_bf > cutofftemp]
otemps_bf = otemps_bf[otemps_bf > cutofftemp]


# Bruker
bruker_name = "/home/karajan/uni/master/ma/analyse/data/crn/data/T2/bruker_t2.data"
bruker_data = np.loadtxt(bruker_name)
bruker_temp = bruker_data[:, 1]
bruker_t2 = bruker_data[:, 3]
bruker_t2_err = bruker_data[:, 4]
bruker_beta = bruker_data[:, 5]
bruker_beta_err = bruker_data[:, 6]



# %%
# gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
# plt.figure(figsize=(6, 6))
f, (pltbeta, pltt2) = plt.subplots(
    2, 1, sharex=True, gridspec_kw={"height_ratios": [1, 4]})
f.subplots_adjust(hspace=0.1)
f.set_size_inches(6, 5)


pltt2.errorbar(
    otemps,
    ot2,
    yerr=ot2err,
    linestyle="None",
    color="tab:blue",
    marker="o",
    label="$f_{L}(\\mathrm{OBI}) = 97.2$ MHz")

pltt2.errorbar(
    otemps_bf,
    ot2_bf,
    yerr=ot2err_bf,
    linestyle="None",
    color="tab:blue",
    marker="*",
    # label="$f_{L}(\\mathrm{OBI}) = 97.2$ MHz"
)

pltt2.errorbar(
    bruker_temp,
    bruker_t2,
    yerr=bruker_t2_err,
    linestyle="None",
    color="tab:red",
    marker="D",
    label="$f_{L}(\\mathrm{Bruker}) = 131.0$ MHz")

pltt2.set_yscale("log")
# pltt2.set_xlim(220, 445)
# pltt2.set_ylim(30, 10**6)
pltt2.tick_params(which="both", direction="in", top=True, right=True)
pltt2.set_xlabel("Temperatur [K]")
pltt2.set_ylabel("$T_2$ [s]")
pltt2.legend(loc=3)


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

# pltbeta.set_xlim(220, 445)
pltbeta.set_ylim(0.5, 2.1)
pltbeta.tick_params(
    which="both",
    direction="in",
    top=True,
    right=True,
    labelbottom=False,
    labeltop=True)
pltbeta.set_ylabel("$\\beta$")

save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/t2")
