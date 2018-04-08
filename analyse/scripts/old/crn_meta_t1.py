import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob, os
from nmr_lib import *


def plot_setup():
    # plt.style.use("ggplot")
    plt.grid(True)
    plt.yscale("log")
    plt.xlabel("Temperatur [K]")
    # plt.ylabel("1/$T_1$ [1/s]")
    plt.ylabel("beta")
    # plt.title("CRN $T_1$ vgl. mit Zürn Daten")
    plt.title("CRN $\\beta_{T_1}$")

plot_setup()


def do_t1():
    fn_extern = "/home/jens/Documents/lit/CRN/zürn paper daten/t1-zuern.data"
    labels_extern = "Zürn T1"
    data = np.loadtxt(fn_extern)
    temp = data[:,0]
    t1 = data[:,1]
    # plot_data(temp, 1/t1, label=labels_extern[i], marker=".")
    plt.scatter(temp, 1/t1, label=labels_extern, marker=".", color="y")


    # t1dir = "/home/jens/Documents/projekte/crn/data/T1"
    # for i, file_name in enumerate(sorted(glob.glob(t1dir + "/*.data"))):
    t1dir = [
            "/home/jens/Documents/projekte/crn/data/T1/t1_230K_280K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_280K_290K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_300K_310K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_310K_355K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_342K_380K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_360K_440K_170807.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_270K_330K.data",
            ]
    for i, file_name in enumerate(t1dir):
        data = np.loadtxt(file_name)
        temp = data[:,1]
        t1 = data[:,3]
        t1_err = data[:,4]

        plt.errorbar(temp, 1/t1, yerr=t1_err/t1**2, fmt='x')

    save_plot("/home/jens/Documents/projekte/crn/data/T1/t1")


def do_beta():
    plt.yscale("linear")
    # t1dir = "/home/jens/Documents/projekte/crn/data/T1"
    # for i, file_name in enumerate(sorted(glob.glob(t1dir + "/*.data"))):
    t1dir = [
            "/home/jens/Documents/projekte/crn/data/T1/t1_230K_280K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_280K_290K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_300K_310K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_310K_355K.data",
            # "/home/jens/Documents/projekte/crn/data/T1/t1_342K_380K.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_360K_440K_170807.data",
            "/home/jens/Documents/projekte/crn/data/T1/t1_270K_330K.data",
            ]
    for i, file_name in enumerate(t1dir):
        data = np.loadtxt(file_name)
        temp = data[:,1]
        beta = data[:,5]
        beta_err = data[:,6]

        beta = beta[temp < 380]
        beta_err = beta_err[temp < 380]
        temp = temp[temp < 380]

        plt.errorbar(temp, beta, yerr=beta_err, fmt='x')

    save_plot("/home/jens/Documents/projekte/crn/data/T1/t1_beta")


# plt.xlim(200, 500)
# plt.ylim(10**1, 10**6)
# save_plot("/home/jens/Documents/projekte/crn/170817/plots/t1_vergleich")

do_t1()
# do_beta()

show_plot()

