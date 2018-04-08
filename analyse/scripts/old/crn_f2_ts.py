import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob, os

plt.style.use("ggplot")
# plt.rc("font", **{"family":"sans-serif", "sans-serif":["Roboto"]})


def t1(x, p1, p2, p3, p4):
    return p1 * (1 - 2*np.exp(-(x/p2)**p3)) + p4


def f2t1(x, p1, p4):
    return p1 * (1 - 2*np.exp(-(x/p2)**p3)) + p4


def t2(x, p1, p2, p3, p4):
    return p1 * (1 - 2*np.exp(-(2*x/p2)**p3)) + p4


p2 = 6.662e-3
p3 = 0.642
def f2(x, p4, p5, p6, p7):
    # p1 * (1 - 2*np.exp(-(x/p2)**p3)) + p4
    return p5 + p4 * np.exp(-(x/p6)**p7) * np.exp(-(x/p2)**p3)


def load_data():
    fn_nmr = glob.glob("*.dat.nmr")[0]
    data_nmr = np.loadtxt(fn_nmr)
    tm = data_nmr[:,0]

    fn_ts = glob.glob("*.ts")
    fn_ts = sorted(fn_ts, key=lambda a: int(a.split("_")[-1].split(".")[0]))
    time = np.loadtxt(fn_ts[0])[:,0]

    real = np.empty((len(fn_ts), np.loadtxt(fn_ts[0]).shape[0]))
    imag = np.empty((len(fn_ts), np.loadtxt(fn_ts[0]).shape[0]))
    for i, fn_tsi in enumerate(fn_ts):
        data_tsi = np.loadtxt(fn_tsi)
        real[i] = data_tsi[:,1]
        imag[i] = data_tsi[:,3]

    return tm, time, real, imag


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


def select_data(skip_start, skip_end, tm, echo_real, echo_imag):
    end_index = tm.size - skip_end
    return (tm[skip_start:end_index], echo_real[skip_start:end_index], 
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


def fit_data(tm, cmplx):
    p = curve_fit(f2, tm, cmplx.real, p0=[80, 80, 1e-3, 0.9])[0]
    fit = f2(np.sort(tm), p[0], p[1], p[2], p[3])
    return p, fit


def plot_data_fit(tm, cmplx, p, fit):
    plt.text(4e-2, -1000, "$T_1$ = {} ms".format(np.round(p[1] * 1000, decimals=3)), fontsize=10)
    plt.text(4e-2, -1500, "$\\beta$ = {}".format(np.round(p[2], decimals=3)), fontsize=10)
    # plt.plot(np.sort(tm), fit)

    plot_data(tm, cmplx)


def plot_data(tm, cmplx):
    plt.xscale("log")

    plt.scatter(tm, cmplx.real, label="real", color="b", marker="x")
    plt.scatter(tm, cmplx.imag, label="imag", color="r", marker="x")
    plt.plot(np.sort(tm), f2(np.sort(tm), p[0], p[1], p[2], p[3]))

    plt.legend()
    plt.title("CRN $T_1$ {}".format(directory))

    # plt.savefig("plots/05_tc_notso_nice.pdf", bbox_inches="tight")
    plt.show()



directory = "/home/jens/Documents/projekte/crn/170706/F2COS/1289_CRN_F2_280K"
os.chdir(directory)

tm, time, real, imag = load_data()
echo_real, echo_imag = pick_peaks(time, real, imag)
# tm, echo_real, echo_imag = select_data(0, 0, tm, echo_real, echo_imag)
# cmplx = phase_shift(phase, echo_real, echo_imag)
phase, cmplx = find_phase(-15, 5, 50, echo_real, echo_imag)
print(phase)
p, fit = fit_data(tm, cmplx)
print(p)
plot_data_fit(tm, cmplx, p, fit)
# plot_data(tm, cmplx)


# home_dir = "/home/jens/Documents/projekte/crn/170706/T1F2"
# for directory in sorted(glob.glob(home_dir + "/*/")):
#     os.chdir(directory)
#     print(directory)

#     tm, time, real, imag = load_data()
#     echo_real, echo_imag = pick_peaks(time, real, imag)
#     # tm, echo_real, echo_imag = select_data(0, 0, tm, echo_real, echo_imag)
#     # cmplx = phase_shift(phase, echo_real, echo_imag)
#     phase, cmplx = find_phase(-10, 5, 50, echo_real, echo_imag)
#     print(phase)
#     p, fit = fit_data(tm, cmplx)
#     print(p[1])
#     # plot_data(tm, cmplx, p, fit)
