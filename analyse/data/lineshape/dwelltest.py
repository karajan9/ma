# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from lmfit import Model, Parameters
from scipy.optimize import curve_fit, leastsq
from scipy.special import lambertw
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(
    os.path.abspath("/home/karajan/uni/master/analyse/scripts"))
from nmr_lib import *


# %%
# CRN taus

def tau_c_4er(T):
    tau_co = 5.1e-14
    D = 3.5
    T_VF = 294
    return tau_co * np.exp((D * T_VF) / (T - T_VF))


def T_VFT_4er(tau_c):
    tau_co = 5.1e-14
    D = 3.5
    T_VF = 294
    return (D * T_VF) / (np.log(tau_c / tau_co)) + T_VF


def T_VFT_aug(tau_c):
    tau_co = 1.15e-14
    D = 4.72
    T_VF = 285
    return (D * T_VF) / (np.log(tau_c / tau_co)) + T_VF


def T_mauro_aug(tau_c):
    tau_co = 4.75e-12
    K = 6.93
    C = 2320
    return (C / (lambertw(C * np.log(tau_c / tau_co) / K))).real


# %%
basedir = "/home/karajan/uni/master/analyse/NMRSimulation/test/"
fn1 = basedir + "dwelltest0.5e-6.fid.fft"
fn2 = basedir + "dwelltest0.5e-6!1.74010323e-05.fid.fft"
fn3 = basedir + "lineshape_lifetime!4.41927079e-05!.fid.fft"
fn4 = basedir + "lineshape_lifetime!3.45255530e-07!.fid.fft"

d1 = np.loadtxt(fn1)
d2 = np.loadtxt(fn2)
d3 = np.loadtxt(fn3)
d4 = np.loadtxt(fn4)

# plt.plot(d1[:,1], d1[:,2])
# plt.plot(d2[:,1], d2[:,2])
plt.plot(d3[:,1], d3[:,2])
plt.plot(d4[:,1], d4[:,2])

plt.xlim(-500000, 500000)
