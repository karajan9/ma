import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, leastsq
import glob, os
from nmr_lib import *


def plot_setup():
    # plt.style.use("ggplot")
    # plt.yscale("log")
    plt.grid(True)
    # plt.xlabel("Temperature")
    plt.xlabel("Frequenz [kHz]")
    # plt.ylabel("FWHM [kHz]")
    plt.ylabel("a.u.")
    plt.title("Spektren")

plot_setup()



dirs = [
        # "/home/jens/Documents/projekte/crn/170731/spektren",
        # "/home/jens/Documents/projekte/crn/170807/SPEK",
        "/home/jens/Documents/projekte/crn/170731/joachimsspek",
        
       ]
labels = [
          # "Messung 1",
          # "Messung 2",
          "Joachims Daten",
         ]
# for i, directory in enumerate(dirs):
#     temp, fwhm = load_FWHMs(directory)
#     print(temp, fwhm)
#     plt.scatter(temp, fwhm/1e3, label=labels[i], marker=".")
#     # plot_setup()


file_name = "/home/jens/Documents/projekte/crn/170807/SPEK/spek_fwhm.data"
file_names = [
              "/home/jens/Documents/projekte/crn/170731/spektren/spek_fwhm.data",
              "/home/jens/Documents/projekte/crn/170807/SPEK/spek_fwhm.data",
              "/home/jens/Documents/projekte/crn/170817/SPEK/spek_fwhm.data",
              "/home/jens/Documents/projekte/crn/170828/SPEK/spek_fwhm.data"
             ]
labels = [
          "Messung 1",
          "Messung 2",
          "360-310",
          "Lehrstuhlv."
         ]
# for i, file_name in enumerate(file_names):
#     data = np.loadtxt(file_name)
#     temp = data[:,1]
#     gamma = data[:,5] * 2 / 1e3
#     gamma_err = data[:,6] * 2 / 1e3

#     plt.errorbar(temp, gamma, yerr=gamma_err, label=labels[i], fmt='.')



directory = "/home/jens/Documents/projekte/crn/170817/SPEK"
# labels = ["380 K", "376 K", "373 K", "370 K", "368 K", "366 K", 
#           "364 K", "362 K", "360 K", "357 K", "354 K", "351 K", 
#           "348 K", "345 K", "351 K", "348 K", "345 K", "342 K"]
labels = ["355 K", "350 K", "345 K", "340 K", "330 K", "325 K", 
          "320 K", "315 K", "310 K", "335 K"]
data = np.loadtxt("/home/jens/Documents/projekte/crn/170817/SPEK/spek_fwhm.data")
for i, fn_spek in enumerate(sorted(glob.glob(directory + "/*/*.spec.nmr"))):
    freq, real, imag = load_spek(fn_spek)
    cmplx = to_cmplx(real, imag)

    phase, cmplx = phase_fit(cmplx, 0)
    print(phase)
    plt.plot(freq/1e3, cmplx.real/cmplx.real.max() - i * 0.75, label=labels[i], color="b")
    plot_setup()
    try:
        p_value, p_error, fit = fit_lorentz(freq, cmplx.real/cmplx.real.max())
        # print(experiment_number, phase, p_value[0], p_error[0], p_value[1], p_error[1], p_value[2], p_error[2])
        # plt.plot(freq/1e3, fit - i * 0.75, color="r")
    except Exception as e:
        print(e)
        print(str(experiment_number) + ": fit failed")
    plt.plot([-250, 250], [- i * 0.75, - i * 0.75], color="y")

    gamma = data[:,5] * 2
    idx = (np.abs(cmplx.real[freq<0]/cmplx.real.max() - 0.5)).argmin()
    print(cmplx.real)
    print(freq[idx])
    plt.plot([freq[idx]/1e3, freq[idx]/1e3 + gamma[i]/1e3], [0.5 - i * 0.75, 0.5 - i * 0.75], color="r")

    
plot_limits(xmin=-250, xmax=250)

# save_plot("/home/jens/Documents/projekte/crn/170828/plots/fwhm fit")
show_plot()
