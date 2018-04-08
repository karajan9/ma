import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, leastsq
import glob, os
from nmr_lib import *


def plot_setup():
    # plt.style.use("ggplot")
    # plt.yscale("log")
    plt.grid(True)
    plt.xlabel("Temperature")
    plt.ylabel("FWHM [kHz]")
    plt.title("Spektrum FWHM")

plot_setup()



file_name = "/home/jens/Documents/projekte/crn/170807/SPEK/spek_fwhm_test.data"

fn_t1 = "/home/jens/Documents/projekte/crn/170807/T1/t1_test.data"

data_t1 = np.loadtxt(fn_t1)
id_t1 = data_t1[:,0]
temp = data_t1[:,1]
t1 = data_t1[:,3]
t1_err = data_t1[:,4]

data = np.loadtxt(file_name)
ids = data[:,0]
temp = data[:,1]
gamma = data[:,5]
gamma_err = data[:,6]


plt.errorbar(temp, gamma*2, yerr=gamma_err, label="gamma", fmt='.')
plt.errorbar(temp[6:], (gamma[6:]-1/t1[6:]/2/np.pi)*2, label="gamma -  1/(t1*2*pi)", fmt='o')
plt.errorbar(temp[6:], 1/t1[6:]/3, label="1/t1/3", fmt='x')


file_name = "/home/jens/Documents/projekte/crn/170817/SPEK/spek_fwhm.data"

fn_t1 = "/home/jens/Documents/projekte/crn/170817/T1/t1.data"

data_t1 = np.loadtxt(fn_t1)
id_t1 = data_t1[:,0]
temp = data_t1[:,1]
t1 = data_t1[:,3]
t1_err = data_t1[:,4]

data = np.loadtxt(file_name)
ids = data[:,0]
temp = data[:,1]
gamma = data[:,5]
gamma_err = data[:,6]


plt.errorbar(temp, gamma*2, yerr=gamma_err, label="gamma", fmt='.')
plt.errorbar(temp, (gamma-1/t1/2/np.pi)*2, yerr=gamma_err, label="gamma -  1/(t1*2*pi)", fmt='o')
plt.errorbar(temp, 1/t1/3, label="1/t1/3", fmt='x')

save_plot("/home/jens/Documents/projekte/crn/170817/plots/spek t1 neu")
show_plot()
