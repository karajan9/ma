# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath("/home/karajan/uni/master/analyse/scripts"))
from nmr_lib import *


# %%
def plot_setup():
    plt.gcf().clear()
    # plt.style.use("ggplot")
    plt.grid(True)
    plt.yscale("log")
    plt.xlabel("Temperatur [K]")


# %%
def do_t2():
    plot_setup()
    plt.ylabel("$T_2$ [s]")
    plt.title("CRN $T_2$")

    # home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T2/"
    # fn = [
    #     home_dir + "t2_230K_280K.data",
    #     home_dir + "t2_300K_310K.data",
    #     home_dir + "t2_360K_440K.data",
    #     home_dir + "t2_342K_380K.data",
    #     home_dir + "t2_270K_330K.data",
    #     home_dir + "t2_305K_345K.data",
    #     home_dir + "t2_305K_325K.data",
    # ]

    # for i, file_name in enumerate(fn):
    #     data = np.loadtxt(file_name)
    #     temp = data[:, 1]
    #     t2 = data[:, 3]
    #     t2_err = data[:, 4]

    #     plt.errorbar(temp, t2, yerr=t2_err, fmt='x')

    home_dir = "/home/karajan/uni/master/analyse/data"
    t1files = glob.glob(home_dir + "/*/T2/T2_*.data")

    for file in t1files:
        data = np.loadtxt(file)
        temp = data[:, 1]
        t2 = data[:, 3]
        t2_err = data[:, 4]

        plt.errorbar(temp, t2, yerr=t2_err, fmt='x')

    save_plot(plt, "/home/karajan/uni/master/analyse/plots/T2/t2_neu")


def do_beta():
    plot_setup()
    plt.yscale("linear")
    plt.xlabel("Temperatur [K]")
    plt.ylabel("$\\beta_{T_2}$ [s]")
    plt.title("CRN $\\beta_{T_2}$")


    # home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T2/"
    # t2dir = [
    #     home_dir + "t2_230K_280K.data",
    #     home_dir + "t2_300K_310K.data",
    #     home_dir + "t2_360K_440K.data",
    #     home_dir + "t2_342K_380K.data",
    #     home_dir + "t2_270K_330K.data",
    #     home_dir + "t2_305K_345K.data",
    #     home_dir + "t2_305K_325K.data",
    # ]
    # for i, file_name in enumerate(t2dir):
    #     data = np.loadtxt(file_name)
    #     temp = data[:, 1]
    #     beta = data[:, 5]
    #     beta_err = data[:, 6]

    #     # beta = beta[temp < 380]
    #     # beta_err = beta_err[temp < 380]
    #     # temp = temp[temp < 380]

    #     plt.errorbar(temp, beta, yerr=beta_err, fmt='x')

    home_dir = "/home/karajan/uni/master/analyse/data"
    t1files = glob.glob(home_dir + "/*/T2/T2_*.data")

    for file in t1files:
        data = np.loadtxt(file)
        temp = data[:, 1]
        beta = data[:, 5]
        beta_err = data[:, 6]

        plt.errorbar(temp, beta, yerr=beta_err, fmt='x')

    save_plot(plt, "/home/karajan/uni/master/analyse/plots/T2/t2_beta_neu")


do_t2()
do_beta()


# %%
def do_t2_bruker():
    plot_setup()
    plt.ylabel("$T_2$ [s]")
    plt.title("CRN $T_2$ Bruker")

    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T2/"
    fn = [
        home_dir + "t2_230K_280K.data",
        home_dir + "t2_300K_310K.data",
        home_dir + "t2_360K_440K.data",
        home_dir + "t2_342K_380K.data",
        home_dir + "t2_270K_330K.data",
        home_dir + "t2_305K_345K.data",
        home_dir + "t2_305K_325K.data",
    ]

    for i, file_name in enumerate(fn):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        t2 = data[:, 3]
        t2_err = data[:, 4]

        plt.errorbar(temp, t2, yerr=t2_err, fmt='.', color="blue")

    bruker_name = home_dir + "bruker_t2.data"
    bruker_data = np.loadtxt(bruker_name)
    bruker_temp = bruker_data[:, 1]
    bruker_t2 = bruker_data[:, 3]
    bruker_t2_err = bruker_data[:, 4]
    plt.errorbar(bruker_temp, bruker_t2, yerr=bruker_t2_err, fmt='s', color="red", label="Bruker")
    plt.legend(loc=3)

    save_plot(plt, "/home/karajan/uni/master/analyse/plots/BRUKER/bruker_t2")


do_t2_bruker()


# %%
def do_t2_beta_bruker():
    plot_setup()
    plt.yscale("linear")
    plt.ylabel("beta")
    plt.title("CRN Bruker $\\beta_{T_2}$")

    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T2/"
    fn = [
        home_dir + "t2_230K_280K.data",
        home_dir + "t2_300K_310K.data",
        home_dir + "t2_360K_440K.data",
        home_dir + "t2_342K_380K.data",
        home_dir + "t2_270K_330K.data",
        home_dir + "t2_305K_345K.data",
        home_dir + "t2_305K_325K.data",
    ]

    for i, file_name in enumerate(fn):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        t2 = data[:, 5]
        t2_err = data[:, 6]

        plt.errorbar(temp, t2, yerr=t2_err, fmt='.', color="blue")

    bruker_name = home_dir + "bruker_t2.data"
    bruker_data = np.loadtxt(bruker_name)
    bruker_temp = bruker_data[:, 1]
    bruker_t2 = bruker_data[:, 5]
    bruker_t2_err = bruker_data[:, 6]
    plt.errorbar(bruker_temp, bruker_t2, yerr=bruker_t2_err, fmt='s', color="red", label="Bruker")
    plt.legend(loc=3)

    save_plot(plt, "/home/karajan/uni/master/analyse/plots/BRUKER/bruker_t2beta")


do_t2_beta_bruker()

# plt.xlim(300, 460)
# show_plot()
