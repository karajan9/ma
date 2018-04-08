import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, leastsq
import glob, os
from nmr_lib import *


def plot_setup():
    # plt.style.use("ggplot")
    plt.grid(True)

plot_setup()


def plot_spek(directory):
    plt.xscale("linear")
    plt.xlabel("Frequenz [Hz]")
    plt.ylabel("Amplitude [a.u.]")
    plt.title("Spektren")

    spek_dirs = sorted(glob.glob(directory + "/*/"))
    for i, spek_dir in enumerate(spek_dirs):
        tau = get_tau(spek_dir)
        freq, real, imag = load_spek(spek_dir)
        real /= np.max(real)

        plt.plot(freq, real, label=tau)


def fwhm(directory):
    plt.xscale("log")
    plt.xlabel("log(tau/s)")
    plt.ylabel("FWHM [kHz]")
    plt.title("FWHM")

    temps, fwhms = load_FWHMs(directory)
    taus = np.zeros(fwhms.shape)

    spek_dirs = sorted(glob.glob(directory + "/*/"))[:fwhms.shape[0]]
    for i, spek_dir in enumerate(spek_dirs):
        taus[i] = get_tau(spek_dir)
        temp = get_temperature(spek_dir)

    p_value, p_error, fit = fit_spek(taus, fwhms, p0=[40e3, 1e-3, 1.0, 1e3])
    print(temp, p_value[1], p_error[1],
                p_value[2], p_error[2],
                p_value[0], p_error[0],
                p_value[3], p_error[3])

    x = np.geomspace(1e-5, 1e-2)
    fit = kohlrausch(x, p_value[0], p_value[1], p_value[2], p_value[3])
    plt.plot(taus, fwhms/1e3, label=temp)
    plt.plot(x, fit/1e3)

    show_plot()


def diff(directory):
    plt.xscale("log")
    plt.xlabel("log(tau/s)")
    plt.ylabel("Difference [nomalized]")
    plt.title("Difference at FWHM")

    spek_dirs = sorted(glob.glob(directory + "/*/"))
    freq, comparison_real, comparison_imag = load_spek(spek_dirs[0])

    temps, fwhms = load_FWHMs(directory)
    taus = np.zeros(fwhms.shape)
    left = np.zeros(fwhms.shape)
    right = np.zeros(fwhms.shape)

    temp = get_temperature(spek_dirs[0])

    for i, spek_dir in enumerate(spek_dirs[:fwhms.shape[0]]):
        taus[i] = get_tau(spek_dir)
        left[i], right[i] = calc_spek_diff(comparison_real, spek_fn=spek_dir)

    plt.plot(taus, left, color="k", ls="--")
    # plt.plot(taus, right, label=temp, ls="-.")
 

def diff2(directory):
    plt.xscale("log")
    plt.xlabel("log(tau/s)")
    plt.ylabel("Difference [nomalized]")
    plt.title("Difference at FWHM")

    spek_dirs = sorted(glob.glob(directory + "/*/"))

    taus, left = np.loadtxt(glob.glob(directory + "/*.data")[0], unpack=True)
    temp = get_temperature(spek_dirs[0])

    # bounds=(-np.inf, np.inf)
    p_value, p_error, fit = fit_spek(taus, left, p0=[-1.0, 1e-2, 1.0, 1.0],
        bounds=(-np.inf, np.inf))
    print(temp, p_value[1], p_error[1],
                p_value[2], p_error[2],
                p_value[0], p_error[0],
                p_value[3], p_error[3])

    # plt.plot(taus, fit)
    plt.plot(taus, left, label=temp)
    # plt.plot(taus, right*2, label=temp, ls="-.")


# plot_limits(xmin=-250, xmax=250)


print("# Messungs_ID Temperature[K] tau[s] tau_error Beta Beta_error M0 M0_error Moff Moff_err")

all_dirs = "/home/jens/Documents/projekte/crn/170918/SPEKkombiniert/temp_abh"
for i, directory in enumerate(sorted(glob.glob(all_dirs + "/*/"))):
    # fwhm(directory)
    diff(directory)
    diff2(directory)
    # pass

# diff("/home/jens/Documents/projekte/crn/170918/SPEKkombiniert/temp_abh/305K")

save_plot("/home/jens/Documents/projekte/crn/170918/plots/diffleft")
show_plot()

