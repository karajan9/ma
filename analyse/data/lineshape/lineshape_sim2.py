# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
from lmfit import Model, Parameters
from scipy.optimize import curve_fit, leastsq
from scipy.special import lambertw

# %%
directory = "/home/karajan/uni/master/analyse/data/lineshape/run2"
fn_fft = glob.glob(directory + "/*.fid.fft")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))

fn_apo = glob.glob(directory + "/*.fid.dat.spec.nmr")
fn_apo = sorted(fn_apo, key=lambda a: float(a.split("!")[1]))

labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))
labels

# %%


def tau_c_4er(T):
    tau_co = 5.1e-14
    D = 3.5
    T_VF = 294
    return tau_co * np.exp((D * T_VF) / (T - T_VF))


def T_VFT_4er(tau_c):
    tau_co = 5.1e-14
    D = 3.5
    T_VF = 294
    return (D * T_VF) / (np.log(tau_c / tau_co)) + T_VF


def T_VFT_aug(tau_c):
    tau_co = 1.15e-14
    D = 4.72
    T_VF = 285
    return (D * T_VF) / (np.log(tau_c / tau_co)) + T_VF


def T_mauro_aug(tau_c):
    tau_co = 4.75e-12
    K = 6.93
    C = 2320
    return C / (lambertw(C * np.log(tau_c / tau_co) / K))


Ts = [
    310.,
    320.,
    325.,
    330.,
    335.,
    340.,
    345.,
    350.,
    355.,
    360.,
    365.,
]
Ts = np.array(Ts)
tau_c(Ts)

T_mauro_aug(labels).real

# %%
plt.rcParams['figure.figsize'] = (12, 8)
for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    plt.plot(data_fft[:, 1], data_fft[:, 2], label=labels[i])

plt.xlim(-1100, 900)
plt.legend()


# %%
def calc_fullwidth_xmax(x, y, height=0.5):
    spek_argmax = np.argmax(y)
    spek_max = y[spek_argmax]
    xmax = spek_max * height

    lower = np.argmin(np.abs(y[:spek_argmax] - xmax))
    higher = np.argmin(np.abs(y[spek_argmax:] - xmax))

    x_diff = x[spek_argmax + higher] - x[lower]
    return x_diff


fwhms = np.empty(11)
print(fwhms)
for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    fwhm = calc_fullwidth_xmax(data_fft[:, 1], data_fft[:, 2], height=0.2)
    print(fwhm, labels[i])
    fwhms[i] = fwhm
fwhms[0] = 30
fwhms[1] = 30

plt.plot(T_mauro_aug(labels).real, fwhms[::-1])


# %%
def plot_fwhm():
    plt.xlabel("Temperatur [K]")
    # plt.ylabel("FWHM [kHz]")
    plt.ylabel("FWHM [kHz]")
    plt.title("FWHM")
    base_dir = "/home/karajan/uni/master/analyse/data/crn/data/SPEK/"
    file_names = [
        base_dir + "spek_230K_280K.data",
        base_dir + "spek_270K_330K.data",
        base_dir + "spek_300K_310K.data",
        base_dir + "spek_310K_355K.data",
        base_dir + "spek_342K_380K.data",
        base_dir + "spek_360K_440K.data",
        base_dir + "spek_305K_345K.data",
    ]

    for i, file_name in enumerate(file_names):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        gamma = data[:, 5] * 2 / 1e3
        gamma_err = data[:, 6] * 2 / 1e3

        plt.errorbar(temp, gamma, yerr=gamma_err, color="tab:blue", fmt='.')


plot_fwhm()
plt.scatter(T_VFT_aug(labels), fwhms / 50.0, color="r")

# 1e-5 oder kürzer: 10 Mio
# 2.6e-4: 100 Mio
# länger: 1 Mrd
