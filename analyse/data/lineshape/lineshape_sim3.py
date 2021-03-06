# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from lmfit import Model, Parameters
from scipy.optimize import curve_fit, leastsq
from scipy.special import lambertw
import sys
home_dir = "/home/karajan/uni/master/ma/analyse"
sys.path.append(os.path.abspath("/home/karajan/uni/master/ma/analyse/scripts"))
from nmr_lib import *



# %%
# CRN taus

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
    return (C / (lambertw(C * np.log(tau_c / tau_co) / K))).real


# %%
Ts = np.array([360., 363., 364., 365., 366., 367., 370., 380.,])
# Ts = np.array([225., 250., 275., 300.,])
tau_c_4er(Ts)

T_VFT_4er(tau_c_4er(Ts)*10)

# labels

# %%
# rename and rescale
directory = "/home/karajan/uni/master/analyse/data/lineshape/rename"
fn_fft = glob.glob(directory + "/*.fid.dat")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    data_fft[:, 1] /= 10000
    print(data_fft)
    np.savetxt(
        fn,
        data_fft,
        delimiter='	',
        header='i	omega	re	im	r	phi',
        footer='',
        comments='#')

# %%
# invert
directory = "/home/karajan/uni/master/analyse/data/lineshape/invert"
fn_fft = glob.glob(directory + "/*.fid.fft")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    data_fft[:, 2] = data_fft[:, 2][::-1]
    data_fft[:, 3] = data_fft[:, 3][::-1] * -1
    np.savetxt(
        fn,
        data_fft,
        delimiter='	',
        header='i	omega	re	im	r	phi',
        footer='',
        comments='#')

# %%
# invert
directory = "/home/karajan/uni/master/ma/analyse/data/lineshape/run11/inv"
fn_fft = glob.glob(directory + "/*.fid.dat")
print(fn_fft)
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    data_fft[:, 1] *= 1
    data_fft[:, 2] = data_fft[:, 2] * -1
    print(data_fft[:, 2])
    # data_fft[:,3] = data_fft[:,3][::-1]
    np.savetxt(
        fn,
        data_fft,
        delimiter='	',
        header='i	omega	re	im	r	phi',
        footer='',
        comments='#')

# %%
# plot fft

# directory = "/home/karajan/uni/master/analyse/data/lineshape/run5 dw05"
directory = "/home/karajan/uni/master/analyse/data/lineshape/joined"
# directory = "/home/karajan/uni/master/analyse/data/lineshape/run6"
# directory = "/home/karajan/uni/master/analyse/data/lineshape/rename"
fn_fft = glob.glob(directory + "/*.spec.nmr")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

# plt.rcParams['figure.figsize'] = (8, 12)
for i, fn in enumerate(fn_fft[::-1]):
    data_fft = np.loadtxt(fn, comments="!")
    # cmplx = to_cmplx(data_fft[:, 2], data_fft[:, 3])
    # phase, cmplx = phase_fit(cmplx, 0)
    # data = cmplx.real
    data = data_fft[:, 1]
    temp = np.round(T_VFT_4er(labels[::-1][i]), 2)
    # plt.plot(data_fft[:, 1], cmplx.real - 0.1*i, label=temp)
    plt.plot(
        data_fft[:, 0] / 1e3,
        data / np.max(data) - 0.1 * i,
        label="{} K ≙ {:.1e} s".format(np.round(temp, 1), labels[::-1][i]),
        linewidth=1.2)

# directory = "/home/karajan/uni/master/analyse/data/lineshape/run5 dw05"
# fn_fft = glob.glob(directory + "/*.fid.fft")
# fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
# labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

# # # plt.rcParams['figure.figsize'] = (12, 8)
# for i, fn in enumerate(fn_fft):
#     data_fft = np.loadtxt(fn)
#     plt.plot(data_fft[:, 1], data_fft[:, 2], label=labels[i])


plt.gcf().set_size_inches(6, 7)
plt.legend(bbox_to_anchor=(1.02, 1.02), loc="upper left")
plt.xlabel("Frequenz [kHz]")
plt.tick_params(
    which="both", direction="in", top=True, right=True, labelleft=False)
plt.xlim(-100, 100)
plt.xlabel("Frequenz [Hz]")
# plt.ylabel("Magnetisierung (a.u.)")
# plt.title("Simulation Lineshape")

save_plot(plt,
          "/home/karajan/uni/master/ma/analyse/plots/SPEK2/sim_lineshape")




# %%
# plot dofft
directory = "/home/karajan/uni/master/analyse/data/lineshape/run3/noapo"
fn_apo = glob.glob(directory + "/*.fid.dat.spec.nmr")
fn_apo = sorted(fn_apo, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_apo))

# plt.rcParams['figure.figsize'] = (12, 8)
for i, fn in enumerate(fn_apo):
    data_fft = np.loadtxt(fn, comments="!")
    # plt.plot(data_fft[:, 0], data_fft[:, 1], label=labels[i])
    plt.plot(
        data_fft[:, 0],
        data_fft[:, 1] / np.max(data_fft[:, 1]),
        label=labels[i])

plt.xlim(-1000, 500)
# plt.xlim(-500000, 500000)
plt.legend()


# %%
# dofft plotten
# plt.rcParams['figure.figsize'] = (12, 8)
directory = "/home/karajan/uni/master/ma/analyse/data/lineshape/run11"
directory = "/home/karajan/uni/master/ma/analyse/data/lineshape/joined2"
fn_apo = glob.glob(directory + "/*.fid.fft")
fn_apo = sorted(fn_apo, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_apo))
print(fn_apo)
for i, fn in enumerate(fn_apo[::-1]):
    data_fft = np.loadtxt(fn, comments="#")
    plt.plot(data_fft[:, 1]/1e3, data_fft[:, 2], label=labels[i])

# plt.xlim(-5000, 5000)
# plt.xlim(-500000, 500000)
plt.legend()


# %%
# T_1 Korrektur präppen
def get_fwhm_temp_t1():
    t1_dir = "/home/karajan/uni/master/analyse/data/"
    t1_fn = [
        t1_dir + "170828/T1/T1_170828.data",
        t1_dir + "170807/T1/T1_170807.data",
        t1_dir + "170731/T1/T1_170731.data",
        t1_dir + "170906/T1/T1_170906.data",
        t1_dir + "170912/T1/T1_170912.data",
        t1_dir + "170918/T1/T1_170918.data",
        t1_dir + "170710/T1/T1_170710.data",
        t1_dir + "170713/T1/T1_170713.data",
        t1_dir + "170706/T1/T1_170706.data",
        t1_dir + "170817/T1/T1_170817.data",
    ]
    t1temps = np.empty(0)
    t1t1 = np.empty(0)
    for fn in t1_fn:
        data = np.loadtxt(fn)
        t1temps = np.hstack([t1temps, data[:, 1]])
        t1t1 = np.hstack([t1t1, data[:, 3]])

    spek_dir = "/home/karajan/uni/master/analyse/data/crn/data/SPEK/"
    spek_fn = [
        spek_dir + "spek_230K_280K.data",
        spek_dir + "spek_270K_330K.data",
        spek_dir + "spek_300K_310K.data",
        spek_dir + "spek_310K_355K.data",
        spek_dir + "spek_342K_380K.data",
        spek_dir + "spek_360K_440K.data",
        spek_dir + "spek_305K_345K.data",
    ]
    spektemps = np.empty(0)
    spekfwhm = np.empty(0)
    for fn in spek_fn:
        data = np.loadtxt(fn)
        spektemps = np.hstack([spektemps, data[:, 1]])
        spekfwhm = np.hstack([spekfwhm, data[:, 5] * 2])

    # beste T1 Temperatur finden
    spekt1 = np.empty(len(spektemps))
    for i, temp in enumerate(spektemps):
        tempindex = np.argmin(np.abs(t1temps - temp))
        spekt1[i] = t1t1[tempindex]

    spekfwhm_t1 = spekfwhm - 1 / spekt1 / np.pi
    return spektemps, spekfwhm, spekt1, spekfwhm_t1


def get_fwhm_temp_t1_2():
    homedir = "/home/karajan/uni/master/ma/analyse/data/crn/data/"

    fwhm_data = np.loadtxt(homedir + "SPEK/spek_342K_380K.data")
    temp_f = fwhm_data[:, 1]
    gamma = fwhm_data[:, 5]

    t1_data = np.loadtxt(homedir + "T1/t1_342K_380K.data")
    t1 = t1_data[:, 3]
    t1err = t1_data[:, 4]

    fwhm_data2 = np.loadtxt(homedir + "SPEK/spek_360K_440K.data")
    temp_f2 = fwhm_data2[:, 1]
    gamma2 = fwhm_data2[:, 5]

    t1_data2 = np.loadtxt(homedir + "T1/t1_360K_440K_170807.data")
    t1_data2 = np.loadtxt(homedir + "T1/T1_170807_betafest.data")
    t12 = t1_data2[:, 3]
    t12err = t1_data2[:, 4]

    # return (np.hstack([temp_f, temp_f2[6:]]),
    #        np.hstack([(gamma - 1 / t1[1::2] / 2 / np.pi) * 2, (gamma2[6:] - 1 / t12[6:] / 2 / np.pi) * 2]),
    #        np.hstack([t1[1::2], t12[6:]]),
    #        np.hstack([t1err[1::2], t12err[6:]]))
    return (np.hstack([temp_f, temp_f2]),
            np.hstack([(gamma - 1 / t1[1::2] / 2 / np.pi) * 2,
                       (gamma2 - 1 / t12 / 2 / np.pi) * 2]),
            np.hstack([t1[1::2], t12]),
            np.hstack([t1err[1::2], t12err]))


# # test ob richtige Temperaturen gefunden wurden
# plt.yscale("log")
# plt.scatter(t1temps, t1t1)
# plt.scatter(spektemps, spekt1)




# %%
# T1 rausrechnen
spektemps, spekfwhm, spekt1, spekfwhm_t1 = get_fwhm_temp_t1()
spektemps_2, spekfwhm_t1_2, t1, t1err = get_fwhm_temp_t1_2()
t1err /= 1e3
plt.scatter(spektemps, spekfwhm / 1e3, label="Halbwertsbreite")
print(spekfwhm_t1_2**2)
plt.errorbar(
    spektemps_2,
    spekfwhm_t1_2 / 1e3,
    yerr=-t1err/t1**2/np.pi,
    marker="o",
    linestyle="None",
    color="tab:orange",
    label="Halbwertsbreite mit $T_1$-Korrektur")

plt.gca().tick_params(direction="in", top="on", right="on")
plt.ylim(-5, 50)
plt.xlabel("Temperatur [K]")
plt.ylabel("Halbwertsbreite [kHz]")
plt.xlim(325, 430)
# plt.title("CRN FWHM mit $T_1$-Korrektur")
# plt.title("CRN FWHM")
plt.text(330, 2, "$f_{L}(\\mathrm{OBI}) = 97.2$ MHz", fontsize=12)
plt.legend(loc=1)

save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/fwhm_t1")



# %%


def get_bruker_info():
    home_bruker = ["/home/karajan/uni/master/ma/analyse/data/bruker/SPEK"]
    fn_spek = glob.glob(home_bruker[0] + "/*/*.spec.nmr")
    temp = np.array([
        320.0, 340.0, 350.0, 350.0, 360.0, 360.0, 360.0, 370.0, 380.0, 390.0,
        375.0, 385.0, 365.0, 355.0, 345.0, 330.0, 335.0, 320.0, 310.0, 340.0
    ])
    dirs = np.array([
        "3", "4", "5", "6", "7", "8", "9", "10", "11", "16", "19", "22", "27",
        "32", "37", "42", "45", "48", "51", "55"
    ])
    temp = temp[np.argsort(dirs)]
    index = np.argsort(fn_spek)
    temps, fwhms, means, maxims = get_multi_info(home_bruker)
    return temp, fwhms[index], means[index], maxims[index]


def get_fwhm_temp_t1():
    bruker_name = "/home/karajan/uni/master/ma/analyse/data/crn/data/T1/bruker_t1.data"
    bruker_data = np.loadtxt(bruker_name)
    bruker_temp = bruker_data[:, 1]
    bruker_t1 = bruker_data[:, 3]
    bruker_t1err = bruker_data[:, 4]

    t1temps = bruker_temp
    t1t1 = bruker_t1
    t1err = bruker_t1err

    temps, fwhms, means, maxims = get_bruker_info()

    # beste T1 Temperatur finden
    spekt1 = np.empty(len(temps))
    spekt1err = np.empty(len(temps))
    for i, temp in enumerate(temps):
        tempindex = np.argmin(np.abs(t1temps - temp))
        spekt1[i] = t1t1[tempindex]
        spekt1err[i] = t1err[tempindex]

    spekfwhm_t1 = fwhms - 1 / spekt1 / np.pi
    return temps, fwhms, spekt1, spekfwhm_t1, spekt1err


spektemps, spekfwhm, spekt1, spekfwhm_t1, spekt1err = get_fwhm_temp_t1()
# spektemps_2, spekfwhm_t1_2, t1, t1err = get_fwhm_temp_t1_2()
spekt1err /= 1e3
plt.scatter(spektemps, spekfwhm / 1e3, label="Halbwertsbreite")
# plt.scatter(spektemps, spekfwhm_t1 / 1e3, label="Halbwertsbreite")
plt.errorbar(
    spektemps,
    spekfwhm_t1 / 1e3,
    yerr=spekt1err / spekt1**2 / np.pi,
    marker="o",
    linestyle="None",
    color="tab:orange",
    label="Halbwertsbreite mit $T_1$-Korrektur")
plt.text(330, 5, "$f_{L}(\\mathrm{Bruker}) = 131.0$ MHz", fontsize=12)

print(spekt1)

plt.gca().tick_params(direction="in", top="on", right="on")
plt.ylim(0, 40)
plt.xlabel("Temperatur [K]")
plt.ylabel("Halbwertsbreite [kHz]")
plt.xlim(325, 400)
plt.legend(loc=1)
save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/fwhm_t1_bruker")


# %%
# Sim FWHM
# plt.rcParams['figure.figsize'] = (12, 8)
def calc_fullwidth_xmax(x, y, height=0.5):
    spek_argmax = np.argmax(y)
    spek_max = y[spek_argmax]
    xmax = spek_max * height

    lower = np.argmin(np.abs(y[:spek_argmax] - xmax))
    higher = np.argmin(np.abs(y[spek_argmax:] - xmax))

    x_diff = x[spek_argmax + higher] - x[lower]
    return x_diff


directory = "/home/karajan/uni/master/analyse/data/lineshape/run5 dw05"
fn_fft = glob.glob(directory + "/*.fid.fft")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

fwhms = np.empty(len(fn_fft))
print(fwhms)
for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    # cmplx = to_cmplx(data_fft[:, 2], data_fft[:, 3])
    # phase, cmplx = phase_fit(cmplx, 0)
    # ydata = cmplx.real
    ydata = data_fft[:, 2]
    fwhm = calc_fullwidth_xmax(data_fft[:, 1], ydata, height=0.5)
    # print(fwhm, labels[i])
    fwhms[i] = fwhm
# fwhms[0] = 30
# fwhms[1] = 30

plt.scatter(T_VFT_4er(labels), fwhms, label="shift=0")
plt.ylim(-3000, 300000)

directory = "/home/karajan/uni/master/analyse/data/lineshape/rename"
fn_fft = glob.glob(directory + "/*.fid.fft")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

fwhms = np.empty(len(fn_fft))
print(fwhms)
for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    # cmplx = to_cmplx(data_fft[:, 2], data_fft[:, 3])
    # phase, cmplx = phase_fit(cmplx, 0)
    # ydata = cmplx.real
    ydata = data_fft[:, 2]
    fwhm = calc_fullwidth_xmax(data_fft[:, 1], ydata, height=0.5)
    # print(fwhm, labels[i])
    fwhms[i] = fwhm
# fwhms[0] = 30
# fwhms[1] = 30

plt.scatter(T_VFT_4er(labels), fwhms, label="shift=0")
# plt.ylim(-3000, 300000)

directory = "/home/karajan/uni/master/analyse/data/lineshape/run3/noapo"

fn_apo = glob.glob(directory + "/*.fid.dat.spec.nmr")
fn_apo = sorted(fn_apo, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_apo))

maxims = np.empty(len(fn_apo))
means = np.empty(len(fn_apo))
fwhm = np.empty(len(fn_apo))
for i, fn in enumerate(fn_apo):
    maxims[i], means[i], fwhm[i], temp = get_dofft_info(fn)

# FWHM Exp / Sim Vergleich
plt.scatter(T_VFT_4er(labels / 128), fwhm * 100 * 2 * np.pi, color="r")



plot_maxima()


# %%
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
            plt.scatter(temps, fwhm / 1e3, color="tab:blue", marker="v")
        elif kind == "maxim":
            plt.scatter(temps, maxims / 1e3, color="tab:blue", marker="v")
        elif kind == "mean":
            plt.scatter(temps, means / 1e3, color="tab:blue", marker="v")


# %%
# Exp / Sim Vergleich

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


def get_multi_info(dirs):
    maxims = np.empty(0)
    means = np.empty(0)
    fwhms = np.empty(0)
    temps = np.empty(0)
    for spek_dir in dirs:
        fn_spek = glob.glob(spek_dir + "/*/*.spec.nmr")
        for i, fn in enumerate(fn_spek):
            maxim, mean, fwhm, temp = get_dofft_info(fn)
            maxims = np.hstack([maxims, maxim])
            means = np.hstack([means, mean])
            fwhms = np.hstack([fwhms, fwhm])
            temps = np.hstack([temps, temp])

    return temps, fwhms, means, maxims


def get_lorentz():
    home_spek = "/home/karajan/uni/master/analyse/data/"
    spek_fns = [
        home_spek + "crn/data/SPEK/spek_230K_280K.data",
        home_spek + "crn/data/SPEK/spek_270K_330K.data",
        home_spek + "crn/data/SPEK/spek_300K_310K.data",
        home_spek + "crn/data/SPEK/spek_305K_345K.data",
        home_spek + "crn/data/SPEK/spek_310K_355K.data",
        home_spek + "crn/data/SPEK/spek_342K_380K.data",
        home_spek + "crn/data/SPEK/spek_360K_440K.data",
    ]
    maxims = np.empty(0)
    means = np.empty(0)
    fwhms = np.empty(0)
    temps = np.empty(0)
    for spek_fn in spek_fns:
        data = np.loadtxt(spek_fn, comments="#")
        maxims = np.hstack([maxims, data[:, 7]])
        means = np.hstack([means, data[:, 7]])
        fwhms = np.hstack([fwhms, data[:, 5]])
        temps = np.hstack([temps, data[:, 1]])

    return temps, fwhms, means, maxims



def get_bruker_info():
    home_bruker = ["/home/karajan/uni/master/ma/analyse/data/bruker/SPEK"]
    fn_spek = glob.glob(home_bruker[0] + "/*/*.spec.nmr")
    temp = np.array([
        320.0, 340.0, 350.0, 350.0, 360.0, 360.0, 360.0, 370.0, 380.0, 390.0,
        375.0, 385.0, 365.0, 355.0, 345.0, 330.0, 335.0, 320.0, 310.0, 340.0
    ])
    dirs = np.array([
        "3", "4", "5", "6", "7", "8", "9", "10", "11", "16", "19", "22", "27",
        "32", "37", "42", "45", "48", "51", "55"
    ])
    temp = temp[np.argsort(dirs)]
    index = np.argsort(fn_spek)
    temps, fwhms, means, maxims = get_multi_info(home_bruker)
    return temp, fwhms[index], means[index], maxims[index]


def get_obi_info():
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
    return get_multi_info(spek_dirs)


def get_sim_info():
    directory = "/home/karajan/uni/master/ma/analyse/data/lineshape/joined"
    fn_apo = glob.glob(directory + "/*.fid.dat.spec.nmr")
    fn_apo = sorted(fn_apo, key=lambda a: float(a.split("!")[1]))
    stemps = np.array(sorted(float(fn.split("!")[1]) for fn in fn_apo))

    smaxims = np.empty(len(fn_apo))
    smeans = np.empty(len(fn_apo))
    sfwhms = np.empty(len(fn_apo))
    for i, fn in enumerate(fn_apo):
        smaxims[i], smeans[i], sfwhms[i], stemp = get_dofft_info(fn)

    return stemps, sfwhms, smeans, smaxims


stemps, sfwhms, smeans, smaxims = get_sim_info()
otemps, ofwhms, omeans, omaxims = get_obi_info()
btemps, bfwhms, bmeans, bmaxims = get_bruker_info()
ltemps, lfwhms, lmeans, lmaxims = get_lorentz()

stemps = T_VFT_4er(stemps)

# plt.scatter(T_VFT_aug(labels / 128), means / 15, color="g", marker="^")
# plt.scatter(T_mauro_aug(labels / 128), means / 15, color="y", marker="^")
# plt.rcParams['figure.figsize'] = (12, 8)

cutofftemp = 360

# FWHM Exp / Sim Vergleich
plt.scatter(
    otemps,
    ofwhms / 1e3,
    color="tab:blue",
    # facecolors="none",
    marker="o",
    label="OBI")
plt.scatter(
    ltemps[ltemps > cutofftemp],
    2*lfwhms[ltemps > cutofftemp] / 1e3,
    color="tab:orange",
    # facecolors="none",
    marker="*")
plt.scatter(
    T_VFT_4er(stemps),
    sfwhms / 1e3,
    color="tab:green",
    # facecolors="none",
    marker="s",
    label="Simulation")
plt.scatter(
    btemps,
    bfwhms / 1e3,
    color="tab:red",
    # s=45,
    # facecolors="none",
    marker="D",
    label="Bruker")

# plt.xlim(340, 370)
# plt.gcf().set_size_inches(7.5, 5)
plt.gca().tick_params(direction="in", top="on", right="on")
# plt.ylim(-25, 5)
# plt.xlim(350, 370)
plt.xlabel("Temperatur [K]")
plt.ylabel("Halbwertsbreite [kHz]")
plt.legend(loc=3)

# save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/fwhm")

# cutofftemp = 360

# # Mean Exp / Sim Vergleich
# plt.scatter(
#     otemps[otemps < cutofftemp],
#     omeans[otemps < cutofftemp] / 1e3,
#     color="tab:blue",
#     # facecolors="none",
#     marker="o",
#     label="OBI")
# plt.scatter(
#     ltemps[ltemps > cutofftemp],
#     lmeans[ltemps > cutofftemp] / 1e3,
#     color="tab:blue",
#     # facecolors="none",
#     marker="*")


# plt.scatter(
#     btemps,
#     bmeans / 1e3,
#     color="tab:red",
#     # s=45,
#     # facecolors="none",
#     marker="D",
#     label="Bruker")
# plt.scatter(
#     stemps[stemps > cutofftemp],
#     smaxims[stemps > cutofftemp] / 1e3,
#     color="tab:green",
#     # facecolors="none",
#     marker="*")
# plt.scatter(
#     stemps[stemps < cutofftemp],
#     smeans[stemps < cutofftemp] / 1e3,
#     color="tab:green",
#     # facecolors="none",
#     marker="s",
#     label="Simulation")


# # plt.xlim(340, 370)
# # plt.gcf().set_size_inches(9, 6)
# # plt.gcf().set_size_inches(7.5, 5)
# plt.gca().tick_params(direction="in", top="on", right="on")
# plt.ylim(-25, 5)
# # plt.xlim(350, 380)
# plt.xlabel("Temperatur [K]")
# plt.ylabel("Schwerpunkt [kHz]")
# plt.legend(loc=3)

# save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/mean4")





# %%
# Single Spektren Vergleich
# home_dir = "/home/karajan/uni/master/analyse/data/"
# simfn = "lineshape/joined/lineshape_lifetime!3395656250000.0!.fid.dat.spec.nmr"
# simfn = "lineshape/joined/lineshape_lifetime!3.162953125e-05!.fid.dat.spec.nmr"
# simfn = "lineshape/joined/lineshape_lifetime!1.30435562e-05!.fid.dat.spec.nmr"
# sim = np.loadtxt(home_dir + simfn, comments="!")
# expfn = "170828/SPEK/1661_CRN_T2_348K/CRN_T2_348K_1661_1.ts.spec.nmr"
# exp2 = np.loadtxt(home_dir + expfn, comments="!")

# plt.plot(sim[:,0], sim[:,1] / np.max(sim[:,1]) * np.max(exp2[:,1]))
# plt.plot(exp2[:,0], exp2[:,1])
# plt.xlim(-2e5, 2e5)
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/SIM/sim_exp_347")

# home_dir = "/home/karajan/uni/master/analyse/data/"
# simfn = "lineshape/joined/lineshape_lifetime!3395656250000.0!.fid.dat.spec.nmr"
# sim = np.loadtxt(home_dir + simfn, comments="!")
# expfn = "lineshape/exp-sel/dummyfolder/CRN_T2_310K_1504_1.ts.spec.nmr"
# exp2 = np.loadtxt(home_dir + expfn, comments="!")

# plt.plot(sim[:, 0], sim[:, 1] / np.max(sim[:, 1]) * np.max(exp2[:, 1]))
# plt.plot(exp2[:, 0], exp2[:, 1])
# plt.xlim(-2e5, 2e5)
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/SIM/sim_exp_310")


home_dir = "/home/karajan/uni/master/ma/analyse/data/"
simfn = "lineshape/joined/lineshape_lifetime!7.844609375e-10!.fid.fft"
sim = np.loadtxt(home_dir + simfn, comments="#")

phase, cmplx = phase_fit(sim[:, 2] + 1j * sim[:, 3], 0)
plt.plot(sim[:, 1], cmplx.real)
plt.xlim(-2e5, 2e5)
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/SIM/sim_exp_310")




# %%
# Exp Lineshape Vergleich
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")

# plt.rcParams['figure.figsize'] = (12, 12)

home_spek = "/home/karajan/uni/master/analyse/data/"
spek_dirs = [
    # home_spek + "170706/SPEK",
    # home_spek + "170713/SPEK",
    # home_spek + "170731/spektren",
    # home_spek + "170807/SPEK",
    # home_spek + "170817/SPEK",
    # home_spek + "170828/SPEK",
    # home_spek + "170906/SPEK/tau15",
    # home_spek + "170912/SPEK/tau_abh/tau0015",
    home_spek + "lineshape/exp-sel",
]

for spek_dir in spek_dirs:
    fn_spek = np.array(glob.glob(spek_dir + "/*/*.spec.nmr"))
    temps = np.empty(len(fn_spek))

    for i, fn in enumerate(fn_spek):
        temps[i] = get_dofft_info(fn)[3]

    tempsind = temps.argsort()
    # print(fn_spek[tempsind])
    fn_spek = fn_spek[tempsind]
    temps = sorted(temps)

    for i, fn in enumerate(fn_spek):
        data_dofft = np.loadtxt(fn, comments="!")
        plt.plot(
            data_dofft[:, 0] / 1e3,
            data_dofft[:, 1] / np.max(data_dofft[:, 1]) - i * 0.2,
            label="{} K".format(np.round(temps[i], 1)),
            linewidth=1.2)


plt.gcf().set_size_inches(6, 5)
# plt.legend(loc='right')
plt.legend(bbox_to_anchor=(1.02, 1.02), loc="upper left")
plt.xlim(-100, 100)
# plt.show()

plt.xlabel("Frequenz [kHz]")
plt.tick_params(
    which="both", direction="in", top=True, right=True, labelleft=False)

# plt.ylabel("Schwerpunkt [kHz]")
save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/SPEK2/spek_lineshape")



# %%
home_bruker = "/home/karajan/uni/master/ma/analyse/data/bruker/SPEK"
fn_spek = glob.glob(home_bruker + "/*/*.spec.nmr")
temp = np.array([
    320.0, 340.0, 350.0, 350.0, 360.0, 360.0, 360.0, 370.0, 380.0, 390.0,
    375.0, 385.0, 365.0, 355.0, 345.0, 330.0, 335.0, 320.0, 310.0, 340.0
])
dirs = np.array([
    "3", "4", "5", "6", "7", "8", "9", "10", "11", "16", "19", "22", "27",
    "32", "37", "42", "45", "48", "51", "55"
])

temp = temp[np.argsort(dirs)]
index = np.argsort(temp)
fn_spek = np.array(sorted(fn_spek))
fn_spek = fn_spek[index]

for i, fn in enumerate(fn_spek):
    data = np.loadtxt(fn, comments="!")
    plt.plot(data[:, 0], data[:, 1] / np.max(data[:, 1]) - 0.2 * i)

plt.xlim(-1e5, 1e5)



# %%
def calcsigma(factor):
    larmorfrequenz = 97.1722e6  # *2*np.pi
    I = 3.0 / 2.0
    hbar = 1.05457180013e-34

    faktorohneeQ = 1 / larmorfrequenz * (3.0 / (2.0 * I * (
        2 * I - 1.0)))**2 * (I * (I + 1) - 3.0 / 4.0)
    sigmasq = factor / faktorohneeQ
    # lasterror = 1.0 / 2.0 * np.sqrt(1 / (faktorohneeQ * factor)) * err
    return np.sqrt(sigmasq) # , lasterror

calcsigma(10000)
