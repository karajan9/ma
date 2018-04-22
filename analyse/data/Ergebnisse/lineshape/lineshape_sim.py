# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
from lmfit import Model, Parameters
from scipy.optimize import curve_fit, leastsq

# %%
directory = "/home/jens/Documents/simulation/Ergebnisse/lineshape/run1"
fn_fft = glob.glob(directory + "/*.fid.fft")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = sorted(float(fn.split("!")[1]) for fn in fn_fft)


# %%
def calc_fullwidth_xmax(x, y, height=0.5):
    spek_argmax = np.argmax(y)
    spek_max = y[spek_argmax]
    xmax = spek_max * height

    lower = np.argmin(np.abs(y[:spek_argmax] - xmax))
    higher = np.argmin(np.abs(y[spek_argmax:] - xmax))

    x_diff = x[spek_argmax + higher] - x[lower]
    return x_diff


# %%
plt.rcParams['figure.figsize'] = (12, 8)
for i, fn in enumerate(fn_fft[0:2]):
    data_fft = np.loadtxt(fn)
    plt.plot(data_fft[:, 1], data_fft[:, 2], label=labels[i])

plt.xlim(-1100, -900)
plt.legend()

# %%
fwhms = np.empty(11)
print(fwhms)
for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    fwhm = calc_fullwidth_xmax(data_fft[:, 1], data_fft[:, 2])
    print(fwhm, labels[i])
    fwhms[i] = fwhm
fwhms[0] = 30
fwhms[1] = 30
# plt.xscale("log")
plt.plot(Ts, fwhms[::-1])

# %%
tau_co = 5.1e-14
D = 3.5
T_VF = 294


def tau_c(T):
    return tau_co * np.exp((D * T_VF) / (T - T_VF))


# %%
def plot_fwhm():
    plt.xlabel("Temperatur [K]")
    # plt.ylabel("FWHM [kHz]")
    plt.ylabel("FWHM [kHz]")
    plt.title("FWHM")
    file_names = [
        "/home/jens/Documents/projekte/crn/data/SPEK/spek_230K_280K.data",
        "/home/jens/Documents/projekte/crn/data/SPEK/spek_270K_330K.data",
        "/home/jens/Documents/projekte/crn/data/SPEK/spek_300K_310K.data",
        "/home/jens/Documents/projekte/crn/data/SPEK/spek_310K_355K.data",
        "/home/jens/Documents/projekte/crn/data/SPEK/spek_342K_380K.data",
        "/home/jens/Documents/projekte/crn/data/SPEK/spek_360K_440K.data",
        "/home/jens/Documents/projekte/crn/data/SPEK/spek_305K_345K.data",
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


plot_fwhm()
plt.scatter(Ts, fwhms[::-1] / 50.0, color="r")



# 1e-5 oder kürzer: 10 Mio
# 2.6e-4: 100 Mio
# länger: 1 Mrd
