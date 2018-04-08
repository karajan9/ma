import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, leastsq
import glob, os
from nmr_lib import *


def plot_setup():
    plt.style.use("ggplot")
    # plt.yscale("log")
    plt.xlabel("Frequenz [kHz]")
    plt.ylabel("u.a.")
    plt.title("Spektrum Fit")
    # plot_limits(xmin=-250, xmax=250)

plot_setup()


def analyze_data(directory, label=""):
    fn_spek = glob.glob(directory + "/*.spec.nmr")[0]
    freq, real, imag = load_spek(fn_spek)
    cmplx = to_cmplx(real, imag)
    phase, cmplx = phase_fit(cmplx, 0)

    temp = get_temperature(directory)
    experiment_number = get_experiment_number(directory)

    plot_setup()
    # scatter_cmplx(reptime, cmplx, label=label)
    plot_data(freq, cmplx.real, label=temp)

    try:
        p_value, p_error, fit = fit_lorentz(freq, cmplx.real)
        print(experiment_number, temp, phase, 
              p_value[0], p_error[0], p_value[1], p_error[1], p_value[2], p_error[2])
        print(p_value)
        plot_fit(freq, fit, p_value)
    except Exception as e:
        print(e)
        print(str(experiment_number) + ": fit failed")

    show_plot()


def quick_analyze(directory):
    fn_spek = glob.glob(directory + "/*.spec.nmr")[0]
    temp, fwhm = load_FWHM(fn_spek)
    experiment_number = get_experiment_number(directory)

    freq, real, imag = load_spek(fn_spek)
    cmplx = to_cmplx(real, imag)
    phase, cmplx = phase_fit(cmplx, 0)
    print(experiment_number, temp, 0.0, 0.0, 0.0, fwhm/2, 0.0, freq[np.argmax(cmplx.real)], 0.0)



print("# Messungs_ID Temperature[K] Phase[degree] a a_error gamma gamma_error x0 x0_error ")

# directory = "/home/jens/Documents/projekte/crn/170731/spektren"
# labels = ["440 K", "435 K", "430 K", "425 K", "420 K", "415 K", 
#           "410 K", "405 K", "400 K", "395 K", "390 K", "385 K", 
#           "380 K", "375 K", "370 K", "365 K", "360 K", "360 K"]
# directory = "/home/jens/Documents/projekte/crn/170817/SPEK"
# labels = ["355 K", "350 K", "345 K", "340 K", "335 K", 
#           "330 K", "325 K", "320 K", "315 K", "310 K"]
directory = "/home/jens/Documents/projekte/crn/170912/SPEK/tau0015"

for i, spek_dir in enumerate(sorted(glob.glob(directory + "/*/"))):
    # analyze_data(spek_dir)
    quick_analyze(spek_dir)

# analyze_data("/home/jens/Documents/projekte/crn/170828/SPEK/1673_CRN_T2_342K", "")


# save_plot("/home/jens/Documents/projekte/crn/170807/plots/spek multi small")
# show_plot()
