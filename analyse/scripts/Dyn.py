# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from lmfit import Model, Parameters
from scipy.optimize import curve_fit, leastsq
from scipy.special import lambertw
import sys
home_dir = "/home/karajan/uni/master/ma/analyse"
sys.path.append(
    os.path.abspath("/home/karajan/uni/master/ma/analyse/scripts"))
from nmr_lib import *



# %%

# T1 OBI
t1_home = home_dir + "/data/crn/data/T1/"
t1dir = [
    t1_home + "t1_230K_280K.data",
    t1_home + "t1_280K_290K.data",
    t1_home + "t1_300K_310K.data",
    t1_home + "t1_310K_355K.data",
    t1_home + "t1_342K_380K.data",
    t1_home + "t1_360K_440K_170807.data",
    t1_home + "t1_270K_330K.data",
    t1_home + "t1_305K_345K.data",
    t1_home + "t1_305K_325K.data",
]
ot1temps = np.empty(0)
ot1 = np.empty(0)
ot1err = np.empty(0)
for file_name in t1dir:
    data = np.loadtxt(file_name)

    ot1temps = np.hstack([ot1temps, data[:, 1]])
    ot1 = np.hstack([ot1, data[:, 3]])
    ot1err = np.hstack([ot1err, data[:, 4]])


# T2 OBI
t2_home = home_dir + "/data/crn/data/T2/"
fn = [
    t2_home + "t2_230K_280K.data",
    t2_home + "t2_300K_310K.data",
    t2_home + "t2_360K_440K.data",
    t2_home + "t2_342K_380K.data",
    t2_home + "t2_270K_330K.data",
    t2_home + "t2_305K_345K.data",
    t2_home + "t2_305K_325K.data",
]
ot2temps = np.empty(0)
ot2 = np.empty(0)
ot2err = np.empty(0)
for file_name in fn:
    data = np.loadtxt(file_name)
    ot2temps = np.hstack([ot2temps, data[:, 1]])
    ot2 = np.hstack([ot2, data[:, 3]])
    ot2err = np.hstack([ot2err, data[:, 4]])


# Zürn beta
beta_fn = home_dir + "/data/crn/data/beta-zuern.data"
betadata = np.loadtxt(beta_fn)
betatemps = betadata[:, 0]
betataus = betadata[:, 1]
fit = np.polyfit(betatemps, betataus, 1)
betafit = np.poly1d(fit)

betatemps = 1000 / betatemps
betataus = 10**betataus


# F2
f2temps = np.array([300., 310., 310.,])
f2sin = np.array([5.26e-3, 3.27e-3, 4.23e-3,])
f2sinerr = np.array([0.94e-3, 1.22e-3, 0.37e-3,])
f2cos = np.array([3.68e-3, 9.95e-3, 3.28e-3,])
f2coserr = np.array([0.63e-3, 4.32e-3, 0.07e-3,])


# Spek Dyn
spek_fn = home_dir + "/data/170918/SPEKkombiniert/temp_abh/fwhm_dynamic.data"
spekdata = np.loadtxt(spek_fn)
spektemps = spekdata[:, 0]
spektaus = spekdata[:, 1]
spektauserr = spekdata[:, 2]



# %%
plt.gcf().set_size_inches(6, 5)
plt.yscale("log")
plt.tick_params(which="both", direction="in", top=True, right=True)

plt.errorbar(
    ot1temps, ot1, yerr=ot1err, linestyle="None", label="$T_1$", marker="x")
plt.errorbar(
    ot2temps,
    ot2,
    yerr=ot2err,
    label="$T_2$",
    linestyle="None",
    marker="x",
    color="tab:purple")

# plt.scatter(betatemps, betataus, label="$\\beta$")
Ts = np.linspace(200, 450, 1000)
plt.plot(Ts, 10**betafit(1000 / Ts), label="$\\beta$-Fit", color="tab:brown")

plt.errorbar(
    f2temps,
    f2sin,
    yerr=f2sinerr,
    label="$F_2$ Sin-Sin",
    linestyle="None",
    marker="s",
    color="tab:green")
plt.errorbar(
    f2temps,
    f2cos,
    yerr=f2coserr,
    label="$F_2$ Cos-Cos",
    linestyle="None",
    color="tab:red",
    marker="s")

plt.errorbar(
    spektemps,
    spektaus,
    yerr=spektauserr,
    label="$t_p$-abh. Spek.",
    linestyle="None",
    marker="D",
    color="tab:orange")


plt.xlim(296, 350)
plt.ylim(5e-5, 2e-2)
plt.xlabel("Temperatur [K]")
plt.ylabel("Zeitkonstante [s]")
plt.legend(bbox_to_anchor=(1.02, 1.02), loc="upper left")


save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/dyn")



# %%
def tau_c_4er(T):
    tau_co = 5.1e-14
    D = 3.5
    T_VF = 294
    return tau_co * np.exp((D * T_VF) / (T - T_VF))


def tau_c_aug(T):
    tau_co = 1.15e-14
    D = 4.72
    T_VF = 285
    return tau_co * np.exp((D * T_VF) / (T - T_VF))


plt.gcf().set_size_inches(4, 3)
T = np.linspace(350, 450, 1000)
omega = np.geomspace(1, 2e12)
# plt.plot(T, tau_c(T), label="Arrhenius")
plt.plot(T, tau_c_4er(T), label="$\\tau_c$ nach [Pim+97]", color="tab:orange")
plt.plot(T, tau_c_aug(T), label="$\\tau_c$ nach [Lun+10]", color="tab:green")
plt.scatter(
    410,
    0.61 / 2 / np.pi / 97.2e6,
    marker='o',
    label="$\\omega \\tau_c = 0.61$",
    color="tab:blue")
plt.legend()
plt.tick_params(which="both", direction="in", top=True, right=True)

plt.yscale("log")
plt.xlabel("Temperatur [K]")
plt.ylabel("$\\tau_c$ [s]")
# plt.title("$\\tau_c$ Vergleich")
# plt.savefig("plots/tau_c_arrhenius_vogel_fulcher.pdf", bbox_inches="tight")

save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/vftau")
