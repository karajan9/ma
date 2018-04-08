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
# home_dir = "/home/karajan/uni/master/analyse/data"
def do_t1():
    plot_setup()
    plt.ylabel("1/$T_1$ [1/s]")
    plt.title("CRN $T_1$ vgl. mit Zürn Daten")

    fn_extern ="/home/karajan/uni/master/analyse/data/crn/data/T1/t1-zuern.data"
    labels_extern = "Zürn T1"
    data = np.loadtxt(fn_extern)
    temp = data[:, 0]
    t1 = data[:, 1]
    # plot_data(temp, 1/t1, label=labels_extern[i], marker=".")
    plt.scatter(temp, 1 / t1, label=labels_extern, marker=".", color="y")

    # t1dir = "/home/jens/Documents/projekte/crn/data/T1"
    # for i, file_name in enumerate(sorted(glob.glob(t1dir + "/*.data"))):

    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T1/"
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

    for i, file_name in enumerate(t1dir):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        t1 = data[:, 3]
        t1_err = data[:, 4]

        plt.errorbar(temp, 1 / t1, yerr=t1_err / t1**2, fmt='x')

    # home_dir = "/home/karajan/uni/master/analyse/data"
    # t1files = glob.glob(home_dir + "/*/T1/T1_*.data")

    # for file in t1files:
    #     data = np.loadtxt(file)
    #     temp = data[:, 1]
    #     t1 = data[:, 3]
    #     t1_err = data[:, 4]

    #     plt.errorbar(temp, 1 / t1, yerr=t1_err / t1**2, fmt='x')

    plt.xlim(200, 450)

    save_plot(plt, "/home/karajan/uni/master/analyse/plots/T1/t1_neu")

# %%
def do_t1_bruker():
    plot_setup()
    plt.ylabel("1/$T_1$ [1/s]")
    plt.title("CRN $T_1$ Bruker")

    fn_extern ="/home/karajan/uni/master/analyse/data/crn/data/T1/t1-zuern.data"
    labels_extern = "Zürn T1"
    data = np.loadtxt(fn_extern)
    temp = data[:, 0]
    t1 = data[:, 1]
    # plot_data(temp, 1/t1, label=labels_extern[i], marker=".")
    plt.scatter(
        temp, 1 / t1, label=labels_extern, marker=".", color="tab:orange")

    # t1dir = "/home/jens/Documents/projekte/crn/data/T1"
    # for i, file_name in enumerate(sorted(glob.glob(t1dir + "/*.data"))):
    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T1/"
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
    for i, file_name in enumerate(t1dir):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        t1 = data[:, 3]
        t1_err = data[:, 4]

        plt.scatter(temp, 1 / t1, color="blue")

    bruker_name="/home/karajan/uni/master/analyse/data/crn/data/T1/bruker_t1.data"
    bruker_data = np.loadtxt(bruker_name)
    bruker_temp = bruker_data[:, 1]
    bruker_t1 = bruker_data[:, 3]
    bruker_t1_err = bruker_data[:, 4]
    plt.errorbar(
        bruker_temp,
        1 / bruker_t1,
        yerr=bruker_t1_err / bruker_t1**2,
        fmt='s',
        color="red",
        label="Bruker")

    plt.xlim(220, 450)
    plt.ylim(30, 10**6)
    plt.legend(loc=2)
    save_plot(plt, "/home/karajan/uni/master/analyse/plots/BRUKER/bruker_t1")


do_t1_bruker()


# %%
def do_bruker_beta():
    plot_setup()
    plt.yscale("linear")
    plt.ylabel("beta")
    plt.title("CRN Bruker $\\beta_{T_1}$")

    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T1/"
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
    for i, file_name in enumerate(t1dir):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        beta = data[:, 5]
        beta_err = data[:, 6]

        beta = beta[temp < 380]
        beta_err = beta_err[temp < 380]
        temp = temp[temp < 380]

        plt.errorbar(temp, beta, yerr=beta_err, fmt='x', color="tab:blue")

    bruker_name="/home/karajan/uni/master/analyse/data/crn/data/T1/bruker_t1.data"
    bruker_data = np.loadtxt(bruker_name)
    bruker_temp = bruker_data[:, 1]
    bruker_beta = bruker_data[:, 5]
    bruker_beta_err = bruker_data[:, 6]
    plt.errorbar(
        bruker_temp,
        bruker_beta,
        yerr=bruker_beta_err,
        fmt='s',
        color="red",
        label="Bruker")
    
    plt.legend(loc=3)
    
    # plt.ylim(0, 1.1)

    # home_dir = "/home/karajan/uni/master/analyse/data"
    # t1files = glob.glob(home_dir + "/*/T1/T1_*.data")

    # for file in t1files:
    #     data = np.loadtxt(file)
    #     temp = data[:, 1]
    #     beta = data[:, 5]
    #     beta_err = data[:, 6]

    #     plt.errorbar(temp, beta, yerr=beta_err, fmt='x')

    save_plot(plt, "/home/karajan/uni/master/analyse/plots/BRUKER/bruker_t1beta")

do_bruker_beta()


# %%
def do_beta():
    plot_setup()
    plt.yscale("linear")
    plt.ylabel("beta")
    plt.title("CRN $\\beta_{T_1}$")

    # t1dir = "/home/jens/Documents/projekte/crn/data/T1"
    # for i, file_name in enumerate(sorted(glob.glob(t1dir + "/*.data"))):
    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/T1/"
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
    for i, file_name in enumerate(t1dir):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        beta = data[:, 5]
        beta_err = data[:, 6]

        beta = beta[temp < 380]
        beta_err = beta_err[temp < 380]
        temp = temp[temp < 380]

        plt.errorbar(temp, beta, yerr=beta_err, fmt='x')
    
    # plt.ylim(0, 1.1)

    # home_dir = "/home/karajan/uni/master/analyse/data"
    # t1files = glob.glob(home_dir + "/*/T1/T1_*.data")

    # for file in t1files:
    #     data = np.loadtxt(file)
    #     temp = data[:, 1]
    #     beta = data[:, 5]
    #     beta_err = data[:, 6]

    #     plt.errorbar(temp, beta, yerr=beta_err, fmt='x')

    save_plot(plt, "/home/karajan/uni/master/analyse/plots/T1/t1_beta")

# do_t1()
# do_beta()


# plt.xlim(200, 500)
# plt.ylim(10**1, 10**6)
# save_plot("/home/jens/Documents/projekte/crn/170817/plots/t1_vergleich")

# %%
# do_t1_bruker()
do_bruker_beta()

# show_plot()
