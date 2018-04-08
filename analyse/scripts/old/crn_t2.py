import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq
import glob, os
from nmr_lib import *


def plot_setup():
    plt.style.use("ggplot")
    plt.xscale("log")
    plt.xlabel("Zeit [s]")
    plt.ylabel("a.u.")
    plt.title("CRN $T_2$")

plot_setup()


def analyze_data(directory, label):
    # directory = "/home/jens/Documents/projekte/crn/170731/T2/1400_CRN_T2_360K"
    os.chdir(directory)

    tau, real, real_err, imag, imag_err = load_data()
    cmplx = to_cmplx(real, imag)
    phase, cmplx = phase_fit(cmplx, 15)
    tau, cmplx = sort_values(tau, cmplx)

    # cmplx, tau = select_data(0, 2e-5, cmplx, indexor=tau, remove=True, select_indexor=True)
    # print(tau)
    cmplx = cmplx[1:]
    tau = tau[1:]
    real_err = real_err[1:]
    # print(tau)

    temp = get_temperature(directory)
    experiment_number = get_experiment_number(directory)

    try:
        p_value, p_error, fit = fit_t2(tau, cmplx, sigma=real_err)
        print(experiment_number, temp, phase, p_value[1], p_error[1],
                                              p_value[2], p_error[2],
                                              p_value[0], p_error[0],
                                              p_value[3], p_error[3])
        plot_fit(tau, fit, p_value)
    except Exception as e:
        print(e)
        print(str(experiment_number) + ": fit failed")

    plot_setup()
    scatter_cmplx(tau, cmplx, label=temp)
    scatter_data(tau, np.abs(cmplx), label=temp)
    show_plot()


    # scatter_cmplx(tau, cmplx, label=label)
    # show_plot()

print("# Messungs_ID Temperature[K] Phase[degree] T2[s] T2_error Beta Beta_error M0 M0_error Moff Moff_err")

# home_dir = "/home/jens/Documents/projekte/crn/170807/T2"
# labels = ["440 K", "435 K", "430 K", "425 K", "420 K", "415 K", 
#           "410 K", "405 K", "400 K", "395 K", "390 K", "385 K", 
#           "380 K", "380 K 2", "375 K", "370 K", "365 K", "360 K"]
home_dir = "/home/jens/Documents/projekte/crn/170918/T2"
# labels = ["355 K", "350 K", "345 K", "340 K", "335 K", 
#           "330 K", "325 K", "320 K", "315 K", "310 K"]

for i, directory in enumerate(sorted(glob.glob(home_dir + "/*/"))):
    analyze_data(directory, "")

# analyze_data("/home/jens/Documents/projekte/crn/170807/T2/1420_CRN_T2_440K", "")

#     text = "{}: $T_1$ = {} ms, $\\beta$ = {}".format(labels[i],
#                                                      np.round(p[1] * 1000, decimals=2),
#                                                      np.round(p[2], decimals=2))
#     plt.text(1.5e-3, -1300 - 200*i, text, fontsize=10)

# show_plot()
# save_plot()
