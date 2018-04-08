import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, leastsq
import glob, os
from nmr_lib import *


def plot_setup():
    # plt.style.use("ggplot")
    # plt.yscale("log")
    plt.grid(True)
    plt.xlabel("Frequenz [kHz]")
    plt.ylabel("u.a.")
    plt.title("Spektrum Vergleich")
    # plot_limits(xmin=-250, xmax=250)

plot_setup()


def analyze_data(directory, label=""):
    fn_spek = glob.glob(directory + "/*/*.spec.nmr")
    labels = ["0.015 ms", "0.1 ms", "0.6 ms"]
    for i, fn in enumerate(sorted(fn_spek)):
        freq, real, imag = load_spek(fn)
        cmplx = to_cmplx(real, imag)
        phase, cmplx = phase_fit(cmplx, 0)
        real = cmplx.real

        plot_setup()
        plot_data(freq, real, label=labels[i])

    plt.xlim(-80000, 80000)
    save_plot("/home/jens/Documents/projekte/crn/170906/plots/spektrum vergleich "+label+" unskaliert")
    show_plot()

    for i, fn in enumerate(sorted(fn_spek)):
        freq, real, imag = load_spek(fn)
        cmplx = to_cmplx(real, imag)
        phase, cmplx = phase_fit(cmplx, 0)
        real = cmplx.real
        real /= np.max(real)

        plot_setup()
        plot_data(freq, real, label=labels[i])
    save_plot("/home/jens/Documents/projekte/crn/170906/plots/spektrum vergleich "+label+" full")
    show_plot()

    for i, fn in enumerate(sorted(fn_spek)):
        freq, real, imag = load_spek(fn)
        cmplx = to_cmplx(real, imag)
        phase, cmplx = phase_fit(cmplx, 0)
        real = cmplx.real
        real /= np.max(real)

        plot_setup()
        plot_data(freq, real, label=labels[i])
        plt.xlim(-80000, 80000)
    save_plot("/home/jens/Documents/projekte/crn/170906/plots/spektrum vergleich "+label)
    show_plot()


    # show_plot()



print("# Messungs_ID Temperature[K] Phase[degree] a a_error gamma gamma_error x0 x0_error ")

directory = "/home/jens/Documents/projekte/crn/170906/SPEK/330K"

# for i, spek_dir in enumerate(sorted(glob.glob(directory + "/*/"))):
#     analyze_data(spek_dir)
    # quick_analyze(spek_dir)

analyze_data(directory, "330K")

# plt.xlim(-80000, 80000)
# show_plot()
