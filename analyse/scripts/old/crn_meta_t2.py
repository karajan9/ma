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
    plt.ylabel("$T_2$ [s]")
    # plt.title("CRN $T_2$")
    plt.title("CRN $\\beta_{T_2}$")

plot_setup()


def do_t2():
    fn = [
          "/home/jens/Documents/projekte/crn/data/T2/t2_230K_280K.data",
          "/home/jens/Documents/projekte/crn/data/T2/t2_300K_310K.data",
          "/home/jens/Documents/projekte/crn/data/T2/t2_360K_440K.data",
          "/home/jens/Documents/projekte/crn/data/T2/t2_342K_380K.data",
          "/home/jens/Documents/projekte/crn/data/T2/t2_270K_330K.data",
         ]

    for i, file_name in enumerate(fn):
        data = np.loadtxt(file_name)
        temp = data[:,1]
        t2 = data[:,3]
        t2_err = data[:,4]

        plt.errorbar(temp, t2, yerr=t2_err, fmt='x')

    save_plot("/home/jens/Documents/projekte/crn/data/T2/t2")



def do_beta():
    plt.yscale("linear")
    # t1dir = "/home/jens/Documents/projekte/crn/data/T1"
    # for i, file_name in enumerate(sorted(glob.glob(t1dir + "/*.data"))):
    t2dir = [
              "/home/jens/Documents/projekte/crn/data/T2/t2_230K_280K.data",
              "/home/jens/Documents/projekte/crn/data/T2/t2_300K_310K.data",
              "/home/jens/Documents/projekte/crn/data/T2/t2_360K_440K.data",
              "/home/jens/Documents/projekte/crn/data/T2/t2_342K_380K.data",
              "/home/jens/Documents/projekte/crn/data/T2/t2_270K_330K.data",
            ]
    for i, file_name in enumerate(t2dir):
        data = np.loadtxt(file_name)
        temp = data[:,1]
        beta = data[:,5]
        beta_err = data[:,6]

        # beta = beta[temp < 380]
        # beta_err = beta_err[temp < 380]
        # temp = temp[temp < 380]

        plt.errorbar(temp, beta, yerr=beta_err, fmt='x')

    save_plot("/home/jens/Documents/projekte/crn/data/T2/t2_beta")



# do_t2()
do_beta()

# plt.xlim(300, 460)
show_plot()

