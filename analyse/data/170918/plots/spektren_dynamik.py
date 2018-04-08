# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, leastsq
import glob, os
import sys
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath("/home/karajan/uni/master/analyse/scripts"))
from nmr_lib import *


# # %%
# def plot_setup():
#     # plt.style.use("ggplot")
#     plt.grid(True)

# plot_setup()


# %%
def plot_spek(directory, take, temp):
    plt.xscale("linear")
    plt.xlabel("Frequenz [Hz]")
    plt.ylabel("Amplitude [normiert]")
    plt.title("Spektren {}K".format(temp))
    plt.xlim(-1e5, 1e5)
    plt.ylim(-0.2, 1.1)

    spek_dirs = np.take(sorted(glob.glob(directory + "/*/")), take)
    for i, spek_dir in enumerate(spek_dirs):
        tau = get_tau(spek_dir)
        freq, real, imag = load_spek(spek_dir)
        real /= np.max(real)

        plt.plot(freq, real, label=tau)
    plt.legend()


# # %%
# def diff(directory):
#     axarr[1].grid(True)
#     # axarr[1].xscale("log")
#     # axarr[1].set_xlabel("log(tau/s)")
#     axarr[1].set_ylabel("Difference [nomalized]")
#     # plt.title("Difference at FWHM")

#     spek_dirs = sorted(glob.glob(directory + "/*/"))
#     freq, comparison_real, comparison_imag = load_spek(spek_dirs[0])

#     temps, fwhms = load_FWHMs(directory)
#     taus = np.zeros(fwhms.shape)
#     left = np.zeros(fwhms.shape)
#     right = np.zeros(fwhms.shape)

#     temp = get_temperature(spek_dirs[0])

#     for i, spek_dir in enumerate(spek_dirs[:fwhms.shape[0]]):
#         taus[i] = get_tau(spek_dir)
#         left[i], right[i] = calc_spek_diff(comparison_real, spek_fn=spek_dir)

#     axarr[1].plot(taus, 0.5 - left)
#     # plt.plot(taus, right, label=temp, ls="-.")
 

# # %%
# def diff2(directory):
#     plt.grid(True)
#     plt.xscale("log")
#     # plt.xlabel("log(tau/s)")
#     plt.ylabel("Difference [nomalized]")
#     # plt.title("Difference at FWHM")

#     spek_dirs = sorted(glob.glob(directory + "/*/"))

#     taus, left = np.loadtxt(glob.glob(directory + "/*.data")[0], unpack=True)
#     temp = get_temperature(spek_dirs[0])

#     # bounds=(-np.inf, np.inf)
#     p_value, p_error, fit = fit_spek(taus, left, p0=[-1.0, 1e-2, 1.0, 0.5],
#         bounds=([-1.5, 0.0, 0.0, 0.4], [-0.5, 1.0, 10.0, 0.6]))
#     print(temp, p_value[1], p_error[1],
#                 p_value[2], p_error[2],
#                 p_value[0], p_error[0],
#                 p_value[3], p_error[3])

#     x = np.geomspace(1e-5, 1e-2)
#     fit = kohlrausch(x, p_value[0], p_value[1], p_value[2], p_value[3])

#     # color = plt.gca()._get_lines.get_next_color()

#     plt.plot(x, fit, color=color)
#     # plt.gca().set_color_cycle(None)
#     plt.scatter(taus, left, label=temp, color=color)
#     # plt.plot(taus, left, label=temp, color=color)
#     # plt.plot(taus, right*2, label=temp, ls="-.")


# %%
# def analyze_data(directory, label):
#     plt.grid(True)
#     plt.xscale("log")
#     plt.xlabel("log(tau/s)")
#     plt.ylabel("Amplitude [a.u.]")
#     plt.title("T_2")
#     # directory = "/home/jens/Documents/projekte/crn/170731/T2/1400_CRN_T2_360K"
#     os.chdir(directory)

#     tau, real, real_err, imag, imag_err = load_data()
#     cmplx = to_cmplx(real, imag)
#     phase, cmplx = phase_fit(cmplx, 15)
#     tau, cmplx = sort_values(tau, cmplx)

#     # cmplx, tau = select_data(0, 2e-5, cmplx, indexor=tau, remove=True, select_indexor=True)
#     # print(tau)
#     cmplx = cmplx[1:]
#     tau = tau[1:]
#     real_err = real_err[1:]
#     # print(tau)

#     temp = get_temperature(directory)
#     experiment_number = get_experiment_number(directory)

#     try:
#         p_value, p_error, fit = fit_t2(tau, cmplx, sigma=real_err)
#         print(experiment_number, temp, phase, p_value[1], p_error[1],
#                                               p_value[2], p_error[2],
#                                               p_value[0], p_error[0],
#                                               p_value[3], p_error[3])
#         # plot_fit(tau, fit, p_value)
#     except Exception as e:
#         print(e)
#         print(str(experiment_number) + ": fit failed")

#     # plot_setup()
#     color = plt.gca()._get_lines.get_next_color()
#     plt.plot(tau, np.abs(cmplx), color=color, label=temp)


# plot_limits(xmin=-250, xmax=250)


# %%
def fwhm(directory):
    plt.grid(True)
    plt.xscale("log")
    plt.xlabel("log($\\tau$/s)")
    # axarr[0].set_ylabel("FWHM [kHz]")
    plt.gcf().set_size_inches(9, 6)
    plt.ylabel("FWHM [kHz]")
    plt.title("FWHM($\\tau$)")

    temps, fwhms = load_FWHMs(directory)
    taus = np.zeros(fwhms.shape)

    spek_dirs = sorted(glob.glob(directory + "/*/"))[:fwhms.shape[0]]
    for i, spek_dir in enumerate(spek_dirs):
        taus[i] = get_tau(spek_dir)
        temp = get_temperature(spek_dir)

    # p_value, p_error, fit = fit_spek(taus, fwhms, p0=[40e3, 1e-3, 1.0, 1e3])

    # print(len(taus))
    newlen = min(len(taus), 11)
    taus = taus[:newlen]
    fwhms = fwhms[:newlen]

    p_value, p_error, fit = fit_spek(taus, fwhms, p0=[40e3, 1e-3, 1.0, 0.0],
        bounds=([10e3, 1e-7, 0.2, 0.0], [60e3, 1e-1, 3.1, 0.1]))

    print(temp, p_value[1], p_error[1],
                p_value[2], p_error[2],
                p_value[0], p_error[0],
                p_value[3], p_error[3])

    color = plt.gca()._get_lines.get_next_color()
    x = np.geomspace(1e-5, 2e-2)
    fit = kohlrausch(x, p_value[0], p_value[1], p_value[2], p_value[3])
    plt.scatter(taus, fwhms/1e3, label=temp, color=color)
    plt.plot(x, fit/1e3, color=color)

    # show_plot()


# %%
def vergleich_fwhm_t2():
    plt.gcf().set_size_inches(9, 6)
    plt.yscale("log")
    plt.xlabel("Temperatur [K]")
    plt.ylabel("log($\\tau$/s)")
    plt.title("$\\tau$ aus FWHM im Vergleich mit $\\tau_{T_2}$")

    t2fn = home_dir + "/data/170912/T2/T2_170912.data"
    t2 = np.loadtxt(t2fn)
    print(t2)
    fwhmfn = home_dir + "/data/170918/SPEKkombiniert/temp_abh/fwhm_dynamic.data"
    fwhm = np.loadtxt(fwhmfn)
    fwhm2fn = home_dir + "/data/170918/SPEKkombiniert/temp_abh/spek_pulslenabh_305K_345K.data"
    fwhm2 = np.loadtxt(fwhm2fn)

    plt.scatter(t2[:,1], t2[:,3], label="$T_2$")
    plt.scatter(fwhm[:,0], fwhm[:,1], label="FWHM")
    # plt.scatter(fwhm2[:,0], fwhm2[:,1], label="fwhm2 tau")

    plt.legend()
    

vergleich_fwhm_t2()
save_plot(plt, "/home/karajan/uni/master/analyse/plots/SPEKDYN/spekdyn_t2")
plt.show()



# %%
plt.gcf().set_size_inches(9, 6)



# %%
# show spectra
dir_305 = home_dir + "/data/170918/SPEKkombiniert/temp_abh/305K"
plot_spek(dir_305, [0, 8, 10, 11], 305)
save_plot(plt, "/home/karajan/uni/master/analyse/plots/SPEKDYN/spekdyn_305K")
plt.show()

dir_325 = home_dir + "/data/170918/SPEKkombiniert/temp_abh/325K"
plot_spek(dir_325, [0, 8, 10, 11], 325)
save_plot(plt, "/home/karajan/uni/master/analyse/plots/SPEKDYN/spekdyn_325K")
plt.show()

dir_345 = home_dir + "/data/170918/SPEKkombiniert/temp_abh/345K"
plot_spek(dir_345, [0, 4, 5, 6], 345)
plt.show()



# %%
# FWHM mit Fits plotten
fwhm_dir = home_dir + "/data/170918/SPEKkombiniert/temp_abh"
print("# Messungs_ID Temperature[K] tau[s] tau_error Beta Beta_error M0 M0_error Moff Moff_err")
for i, directory in enumerate(sorted(glob.glob(fwhm_dir + "/*/"))):
    fwhm(directory)
    
plt.legend()
save_plot(plt, "/home/karajan/uni/master/analyse/plots/SPEKDYN/spekdyn_fits")

# # %%
# plt.gca().set_prop_cycle(None)
# # plt.subplot(312)
# for i, directory in enumerate(sorted(glob.glob(all_dirs + "/*/"))):
#     # plot_spek(directory)
#     diff(directory)

# # %%
# plt.gca().set_prop_cycle(None)
# # plt.subplot(313)
# t2_dir = home_dir + "/data/170912/T2"
# for i, directory in enumerate(sorted(glob.glob(t2_dir + "/*/"))):
#     analyze_data(directory, "")

# plot_spek("/home/jens/Documents/projekte/crn/170918/SPEKkombiniert/temp_abh/305K")

# # %%
# plt.xlabel("log(tau/s)")
# plt.xlim(0, 0.75e-2)
# # plt.ylim(-0.05, 0.7)
# # save_plot("/home/jens/Documents/projekte/crn/170918/plots/spek3")
# show_plot()

