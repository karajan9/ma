import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob, os

plt.style.use("ggplot")
# plt.rc("font", **{"family":"sans-serif", "sans-serif":["Roboto"]})


def t1(x, p1, p2, p3, p4):
    return p1 * (1 - 2*np.exp(-(x/p2)**p3)) + p4


def t2(x, p1, p2, p3, p4):
    return p1 * (1 - 2*np.exp(-(2*x/p2)**p3)) + p4


def load_data():
    fn_nmr = glob.glob("*.dat.nmr")[0]
    data_nmr = np.loadtxt(fn_nmr)
    reptime = data_nmr[:,0]

    fn_ts = glob.glob("*.ts")
    fn_ts = sorted(fn_ts, key=lambda a: int(a.split("_")[-1].split(".")[0]))
    time = np.loadtxt(fn_ts[0])[:,0]

    real = np.empty((len(fn_ts), np.loadtxt(fn_ts[0]).shape[0]))
    imag = np.empty((len(fn_ts), np.loadtxt(fn_ts[0]).shape[0]))
    for i, fn_tsi in enumerate(fn_ts):
        data_tsi = np.loadtxt(fn_tsi)
        real[i] = data_tsi[:,1]
        imag[i] = data_tsi[:,3]

    return reptime, time, real, imag


def pick_peaks(time, real, imag):
    echo_time = 18e-6
    echo_width = 4e-6
    echo_begin = echo_time - echo_width / 2
    echo_end = echo_time + echo_width / 2
    echo_begin_index = np.argmin(np.abs(time - echo_begin))
    echo_end_index = np.argmin(np.abs(time - echo_end))

    echo_real = np.mean(real[:,echo_begin_index:echo_end_index], axis=1)
    echo_imag = np.mean(imag[:,echo_begin_index:echo_end_index], axis=1)
    
    return echo_real, echo_imag


def select_data(skip_start, skip_end, reptime, echo_real, echo_imag):
    end_index = reptime.size - skip_end
    return (reptime[skip_start:end_index], echo_real[skip_start:end_index], 
                                           echo_imag[skip_start:end_index])


def phase_shift(phase, echo_real, echo_imag):
    cmplx = np.empty(echo_real.shape, dtype=complex)
    cmplx.real = echo_real
    cmplx.imag = echo_imag
    cmplx *= np.exp(-2 * np.pi * phase/360 * 1j)

    return cmplx


def find_phase(start, end, steps, echo_real, echo_imag):
    min_imag = np.inf
    min_phase = 0
    min_cmplx = np.empty(echo_real.shape, dtype=complex)
    for phase in np.linspace(start, end, steps):
        cmplx = phase_shift(phase, echo_real, echo_imag)
        img_norm = np.linalg.norm(cmplx.imag)
        if img_norm < min_imag:
            min_imag = img_norm
            min_phase = phase
            min_cmplx = cmplx
    return min_phase, min_cmplx


def fit_data(reptime, cmplx):
    p = curve_fit(t1, reptime, cmplx.real, p0=[80, 1e-3, 0.9, 80])[0]
    fit = t1(np.sort(reptime), p[0], p[1], p[2], p[3])
    return p, fit


# def plot_ts(tm, cmplx):
#     plt.scatter(tm, cmplx.real, label="real", color="b", marker="x")
#     plt.scatter(tm, cmplx.imag, label="imag", color="r", marker="x")
#     # plt.plot(np.sort(tm), f2(np.sort(tm), p[0], p[1], p[2], p[3]))



def plot_fit(reptime, fit, p):
    plt.xscale("log")
    # text = "$T_1$ = {} ms, $\\beta$ = {}".format(np.round(p[1] * 1000, decimals=3), np.round(p[2], decimals=3))
    # plt.text(1e-2, -1000, text, fontsize=10)
    # plt.text(4e-2, -1500, "$\\beta$ = {}".format(np.round(p[2], decimals=3)), fontsize=10)
    plt.plot(np.sort(reptime), fit)


def plot_data(reptime, cmplx, label=""):
    plt.xscale("log")
    plt.scatter(reptime, cmplx.real, label=label, marker="x")
    # plt.scatter(reptime, cmplx.imag, label="imag", color="r", marker="x")


def show_plot():
    plot_descriptions()
    plt.show()


def plot_descriptions():
    plt.legend()
    plt.title("CRN $T_1$, $t_m$ = 1 ms")
    plt.xlabel("Zeit [s]")
    plt.ylabel("a.u.")


def save_plot(fn=None):
    plot_descriptions()
    file_types = ["pdf", "svg"]
    fn_save = "/home/jens/Documents/projekte/crn/170713/plots/t1_300_310"
    if fn is not None:
        fn_save = fn
    for file_type in file_types:
        plt.savefig(fn_save + "." + file_type, bbox_inches="tight")



directory = "/home/jens/Documents/projekte/crn/170731/T1/1399_CRN_T1_360K"
os.chdir(directory)

reptime, time, real, imag = load_data()
echo_real, echo_imag = pick_peaks(time, real, imag)
# reptime, echo_real, echo_imag = select_data(0, 0, reptime, echo_real, echo_imag)
# cmplx = phase_shift(phase, echo_real, echo_imag)
phase, cmplx = find_phase(-15, 5, 50, echo_real, echo_imag)
# print(phase)
# p, fit = fit_data(reptime, cmplx)
# print(p)
# plot_fit(reptime, fit, p)
plot_data(reptime, cmplx)


# home_dir = "/home/jens/Documents/projekte/crn/170713/T1"
# labels = ["300 K, vorher  ",
#           "300 K, nachher",
#           "310 K, vorher  ",
#           "310 K, nachher"]
# for i, directory in enumerate(sorted(glob.glob(home_dir + "/*/"))):
#     os.chdir(directory)
#     print(directory)

#     reptime, time, real, imag = load_data()
#     echo_real, echo_imag = pick_peaks(time, real, imag)
#     # reptime, echo_real, echo_imag = select_data(0, 0, reptime, echo_real, echo_imag)
#     # cmplx = phase_shift(phase, echo_real, echo_imag)
#     phase, cmplx = find_phase(-30, -10, 50, echo_real, echo_imag)
#     # print(phase)
#     p, fit = fit_data(reptime, cmplx)
#     print(phase, p[1], p[2])
#     plot_fit(reptime, fit, p)
#     plot_data(reptime, cmplx, label=labels[i])

#     text = "{}: $T_1$ = {} ms, $\\beta$ = {}".format(labels[i],
#                                                      np.round(p[1] * 1000, decimals=2),
#                                                      np.round(p[2], decimals=2))
#     plt.text(1.5e-3, -1300 - 200*i, text, fontsize=10)

show_plot()
# save_plot()
