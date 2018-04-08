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
    plt.ylabel("1/$T_1$ [1/s]")
    plt.title("CRN $T_1$ vgl. mit Zürn Daten")

plot_setup()

k_b = 1.38e-23

def activation_energy(temp):
    r = 1000000*np.exp(- 67/8.314 * temp)
    print(r)
    return r


x = np.linspace(0, 0.06)
# plot_fit(x, activation_energy(x), 0, label="1")
# plot_fit(x, x * 67, 0, label="2")
# plot_fit(x, np.full(x.shape, 1/67000), 0, label="3")
# plot_fit(x, 1/ x / 67, 0, label="4")

# plt.scatter(3.428, 10**-3.58)
# plt.scatter(3.008, 10**-5.3)

fn_extern = "/home/jens/Documents/lit/CRN/zürn paper daten/t1-zuern.data"
labels_extern = "Zürn T1"
data = np.loadtxt(fn_extern)
temp = data[:,0]
t1 = data[:,1]
# plot_data(temp, 1/t1, label=labels_extern[i], marker=".")
plt.scatter(1/temp, t1, label=labels_extern, marker=".")


fn = "/home/jens/Documents/lit/CRN/zürn paper daten/beta-zuern.data"
labels = "Beta Zürn"
data_b = np.loadtxt(fn)
temp_b = data_b[:,0]
beta_b = data_b[:,1]
plt.scatter(temp_b/1000, 10**beta_b, label=labels, marker="x")

try:
    slope, intercept, std_err = fit_linear(temp_b/1000, np.log(10**beta_b))
    print(slope, intercept, std_err)
    plt.plot(x, 10**(intercept + slope * x))
except Exception as e:
    print(e)
    print("fit failed")





plt.xlim(0, 0.06)
plt.ylim(10**-6, 10**3)
# save_plot("/home/jens/Documents/projekte/crn/170817/plots/t1_vergleich")
show_plot()

