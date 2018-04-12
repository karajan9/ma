# %%
import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit, minimize, leastsq
import glob
# import os
from nmr_lib import load_spek, to_cmplx, phase_fit, fit_lorentz, show_plot,\
                    save_plot


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

# dirs = [
#         # "/home/jens/Documents/projekte/crn/170731/spektren",
#         # "/home/jens/Documents/projekte/crn/170807/SPEK",
#         "/home/jens/Documents/projekte/crn/170731/joachimsspek",

#        ]
# labels = [
#           # "Messung 1",
#           # "Messung 2",
#           "Joachims Daten",
#          ]
# for i, directory in enumerate(dirs):
#     temp, fwhm = load_FWHMs(directory)
#     print(temp, fwhm)
#     plt.scatter(temp, fwhm/1e3, label=labels[i], marker=".")
#     # plot_setup()


def plot_maxima():
    plt.xlabel("Temperatur [K]")
    # plt.ylabel("FWHM [kHz]")
    plt.ylabel("Maximum [kHz]")
    plt.title("Spektren")
    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/SPEK/"
    file_names = [
        home_dir + "spek_230K_280K.data",
        home_dir + "spek_270K_330K.data",
        home_dir + "spek_300K_310K.data",
        home_dir + "spek_310K_355K.data",
        home_dir + "spek_342K_380K.data",
        home_dir + "spek_360K_440K.data",
        home_dir + "spek_305K_345K.data",
    ]
    labels = [
        "Messung 1",
        "Lehrstuhlv 2",
        "F2 Suche",
        " ",
        "Lehrstuhlv 1",
        "Hochtemperatur",
        "Pulslängenabh. Spektren",
    ]
    for i, file_name in enumerate(file_names):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        x0 = data[:, 7] / 1e3
        x0_err = data[:, 8] / 1e3

        plt.errorbar(temp, x0, yerr=x0_err, label=labels[i], fmt='.')


# %%
def get_dofft_info(fn):
    maxim = 0.0
    mean = 0.0
    fwhm = 0.0
    temp = 0.0
    with open(fn, "r") as info_file:
        for line in info_file:
            if "!frequency peak:" in line:
                maxim = float(line.split(" ")[2])
            elif "!mean frequency:" in line:
                mean = float(line.split(" ")[2])
            elif "!FWHM:" in line:
                fwhm = float(line.split(" ")[-2])
            elif "!Sample Temperature" in line:
                temp = float(line.split(" ")[-2])
    return maxim, mean, fwhm, temp


def plot_dofft(kind):
    home_spek = "/home/karajan/uni/master/analyse/data/"
    spek_dirs = [
        home_spek + "170706/SPEK",
        home_spek + "170713/SPEK",
        home_spek + "170731/spektren",
        home_spek + "170807/SPEK",
        home_spek + "170817/SPEK",
        home_spek + "170828/SPEK",
        home_spek + "170906/SPEK/tau15",
        home_spek + "170912/SPEK/tau_abh/tau0015",
    ]
    for spek_dir in spek_dirs:
        fn_spek = glob.glob(spek_dir + "/*/*.spec.nmr")
        print(spek_dir)
        maxims = np.empty(len(fn_spek))
        means = np.empty(len(fn_spek))
        fwhm = np.empty(len(fn_spek))
        temps = np.empty(len(fn_spek))
        for i, fn in enumerate(fn_spek):
            maxims[i], means[i], fwhm[i], temps[i] = get_dofft_info(fn)
        if kind == "fwhm":
            plt.scatter(temps, fwhm / 1e3, color="b", marker="o")
        elif kind == "maxim":
            plt.scatter(temps, maxims / 1e3, color="b", marker="o")
        elif kind == "mean":
            plt.scatter(temps, means / 1e3, color="b", marker="o")


def plot_maxima_bruker():
    plt.xlabel("Temperatur [K]")
    plt.ylabel("Schwerpunkt [kHz]")
    plt.title("Schwerpunkt Bruker")

    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/SPEK/"
    bruker_name = home_dir + "bruker_spek_quick.data"
    bruker_data = np.loadtxt(bruker_name)
    bruker_temp = bruker_data[:, 1]
    bruker_x0 = bruker_data[:, 7] / 1e3
    bruker_x0_err = bruker_data[:, 8] / 1e3
    plt.errorbar(
        bruker_temp,
        bruker_x0,
        yerr=bruker_x0_err,
        color="tab:red",
        label="Bruker",
        fmt='s')


plot_dofft("mean")
plot_maxima_bruker()
plt.legend()
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/BRUKER/bruker_mean")


# %%
def plot_fwhm():
    plt.xlabel("Temperatur [K]")
    # plt.ylabel("FWHM [kHz]")
    plt.ylabel("FWHM [kHz]")
    plt.title("FWHM")
    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/SPEK/"
    file_names = [
        home_dir + "spek_230K_280K.data",
        home_dir + "spek_270K_330K.data",
        home_dir + "spek_300K_310K.data",
        home_dir + "spek_310K_355K.data",
        home_dir + "spek_342K_380K.data",
        home_dir + "spek_360K_440K.data",
        home_dir + "spek_305K_345K.data",
    ]
    labels = [
        "Messung 1",
        "Lehrstuhlv 2",
        "F2 Suche",
        " ",
        "Lehrstuhlv 1",
        "Hochtemperatur",
        "Pulslängenabh. Spektren",
    ]
    for i, file_name in enumerate(file_names):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        gamma = data[:, 5] * 2 / 1e3
        gamma_err = data[:, 6] * 2 / 1e3

        plt.errorbar(
            temp,
            gamma,
            yerr=gamma_err,
            color="tab:blue",
            label=labels[i],
            fmt='.')


# %%
def plot_fwhm_bruker():
    plt.xlabel("Temperatur [K]")
    # plt.ylabel("FWHM [kHz]")
    plt.ylabel("FWHM [kHz]")
    plt.title("Bruker FWHM")
    home_dir = "/home/karajan/uni/master/analyse/data/crn/data/SPEK/"
    file_names = [
        home_dir + "spek_230K_280K.data",
        home_dir + "spek_270K_330K.data",
        home_dir + "spek_300K_310K.data",
        home_dir + "spek_310K_355K.data",
        home_dir + "spek_342K_380K.data",
        home_dir + "spek_360K_440K.data",
        home_dir + "spek_305K_345K.data",
    ]
    labels = [
        "Messung 1",
        "Lehrstuhlv 2",
        "F2 Suche",
        " ",
        "Lehrstuhlv 1",
        "Hochtemperatur",
        "Pulslängenabh. Spektren",
    ]
    for i, file_name in enumerate(file_names):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        gamma = data[:, 5] * 2 / 1e3
        gamma_err = data[:, 6] * 2 / 1e3

        plt.errorbar(temp, gamma, yerr=gamma_err, color="tab:blue", fmt='o')

    bruker_name = home_dir + "bruker_spek_quick.data"
    bruker_data = np.loadtxt(bruker_name)
    bruker_temp = bruker_data[:, 1]
    bruker_gamma = bruker_data[:, 5] * 2 / 1e3
    bruker_gamma_err = bruker_data[:, 6] * 2 / 1e3
    plt.errorbar(
        bruker_temp,
        bruker_gamma,
        yerr=bruker_gamma_err,
        color="tab:red",
        label="Bruker",
        fmt='s')

plt.xlim(370, 400)
plt.ylim(5, 25)
plot_fwhm_bruker()
plt.legend(loc=3)
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/BRUKER/bruker_fwhm")


# %%
def plot_spek_multi():
    directory = "/home/karajan/uni/master/analyse/data/bruker/SPEK"
    # labels = ["380 K", "376 K", "373 K", "370 K", "368 K", "366 K",
    #           "364 K", "362 K", "360 K", "357 K", "354 K", "351 K",
    #           "348 K", "345 K", "351 K", "348 K", "345 K", "342 K"]
    # labels = [
    #     "355 K", "350 K", "345 K", "340 K", "330 K", "325 K", "320 K", "315 K",
    #     "310 K", "335 K"
    # ]
    temp = np.array([
        320.0, 340.0, 350.0, 350.0, 360.0, 360.0, 360.0, 370.0, 380.0, 390.0,
        375.0, 385.0, 365.0, 355.0, 345.0, 330.0, 335.0, 320.0, 310.0, 340.0
    ])
    dirs = np.array([
        "3", "4", "5", "6", "7", "8", "9", "10", "11", "16", "19", "22", "27",
        "32", "37", "42", "45", "48", "51", "55"
    ])
    tempindex = np.argsort(temp)
    temp = temp[tempindex]
    dirs = dirs[tempindex]
    # data = np.loadtxt(
    #     "/home/jens/Documents/projekte/crn/170817/SPEK/spek_fwhm.data")
    for i, fn_spek in enumerate(dirs):
        freq, real, imag = load_spek(
            directory + "/" + fn_spek + "/fid.ts.spec.nmr")
        cmplx = to_cmplx(real, imag)

        phase, cmplx = phase_fit(cmplx, 0)
        print(phase)
        plt.plot(
            freq / 1e3,
            cmplx.real / cmplx.real.max() - i * 0.1,
            label=temp[i])
        plot_setup()
    
    plt.title("Bruker Spektren")
        # try:
        #     p_value, p_error, fit = fit_lorentz(freq,
        #                                         cmplx.real / cmplx.real.max())
        #     # print(experiment_number, phase, p_value[0], p_error[0], p_value[1], p_error[1], p_value[2], p_error[2])
        #     # plt.plot(freq/1e3, fit - i * 0.75, color="r")
        # except Exception as e:
        #     print(e)
        #     print(str(experiment_number) + ": fit failed")
        # plt.plot([-250, 250], [-i * 0.75, -i * 0.75], color="y")

        # gamma = data[:, 5] * 2
        # idx = (np.abs(cmplx.real[freq < 0] / cmplx.real.max() - 0.5)).argmin()
        # print(cmplx.real)
        # print(freq[idx])
        # plt.plot(
        #     [freq[idx] / 1e3, freq[idx] / 1e3 + gamma[i] / 1e3],
        #     [0.5 - i * 0.75, 0.5 - i * 0.75],
        #     color="r")


# plot_limits(xmin=-250, xmax=250)

# plot_maxima()
# save_plot("/home/jens/Documents/projekte/crn/data/SPEK/maxima")
# plot_maxima_bruker()
# save_plot("/home/jens/Documents/projekte/crn/data/SPEK/maxima_bruker")
# plot_fwhm_bruker()
# save_plot("/home/jens/Documents/projekte/crn/data/SPEK/fwhm_bruker")
plt.xlim(-100, 100)
plot_spek_multi()
plt.legend(loc=5)
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/BRUKER/bruker_lineshape")

show_plot()
