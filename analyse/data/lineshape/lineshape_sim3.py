# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from lmfit import Model, Parameters
from scipy.optimize import curve_fit, leastsq
from scipy.special import lambertw
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath("/home/karajan/uni/master/analyse/scripts"))
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
neuelt = np.array([
    336.,
    337.,
    338.,
    339.,
    341.,
    342.,
    343.,
    344.,
])
neuelt2 = np.array([
    334.0,
    335.5,
    336.5,
    337.5,
    338.5,
    342.5,
    344.5,
])

tau_c_4er(347)

# labels


# %%
# rename and rescale
directory = "/home/karajan/uni/master/analyse/data/lineshape/rename"
fn_fft = glob.glob(directory + "/*.fid.dat")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    data_fft[:,1] /= 10000
    print(data_fft)
    np.savetxt(fn, data_fft, delimiter='	', header='i	omega	re	im	r	phi', footer='', comments='#')


# %%
# invert
directory = "/home/karajan/uni/master/analyse/data/lineshape/invert"
fn_fft = glob.glob(directory + "/*.fid.fft")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    data_fft[:,2] = data_fft[:,2][::-1]
    data_fft[:,3] = data_fft[:,3][::-1] * -1
    np.savetxt(fn, data_fft, delimiter='	', header='i	omega	re	im	r	phi', footer='', comments='#')


# %%
# invert
directory = "/home/karajan/uni/master/analyse/data/lineshape/run7"
fn_fft = glob.glob(directory + "/*.fid.dat")
print(fn_fft)
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

for i, fn in enumerate(fn_fft):
    data_fft = np.loadtxt(fn)
    data_fft[:,2] = data_fft[:,2] * -1
    print(data_fft[:,2])
    # data_fft[:,3] = data_fft[:,3][::-1]
    np.savetxt(fn, data_fft, delimiter='	', header='i	omega	re	im	r	phi', footer='', comments='#')


# %%
# plot fft

# directory = "/home/karajan/uni/master/analyse/data/lineshape/run5 dw05"
directory = "/home/karajan/uni/master/analyse/data/lineshape/joined"
# directory = "/home/karajan/uni/master/analyse/data/lineshape/run6"
# directory = "/home/karajan/uni/master/analyse/data/lineshape/rename"
fn_fft = glob.glob(directory + "/*.spec.nmr")
fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

plt.rcParams['figure.figsize'] = (8, 12)
for i, fn in enumerate(fn_fft[::-1]):
    data_fft = np.loadtxt(fn, comments="!")
    # cmplx = to_cmplx(data_fft[:, 2], data_fft[:, 3])
    # phase, cmplx = phase_fit(cmplx, 0)
    # data = cmplx.real
    data = data_fft[:, 1]
    temp = np.round(T_VFT_4er(labels[::-1][i]), 2)
    # plt.plot(data_fft[:, 1], cmplx.real - 0.1*i, label=temp)
    plt.plot(data_fft[:, 0], data / np.max(data) - 0.1*i, label=temp)


# directory = "/home/karajan/uni/master/analyse/data/lineshape/run5 dw05"
# fn_fft = glob.glob(directory + "/*.fid.fft")
# fn_fft = sorted(fn_fft, key=lambda a: float(a.split("!")[1]))
# labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_fft))

# # plt.rcParams['figure.figsize'] = (12, 8)
# for i, fn in enumerate(fn_fft):
#     data_fft = np.loadtxt(fn)
#     plt.plot(data_fft[:, 1], data_fft[:, 2], label=labels[i])


plt.xlim(-75000, 50000)
plt.legend(loc=5)

plt.xlabel("Frequenz [Hz]")
# plt.ylabel("Magnetisierung (a.u.)")
plt.title("Simulation Lineshape")

save_plot(plt, "/home/karajan/uni/master/analyse/plots/SIM/sim_lineshape2")


# %%
# plot dofft
directory = "/home/karajan/uni/master/analyse/data/lineshape/run3/noapo"
fn_apo = glob.glob(directory + "/*.fid.dat.spec.nmr")
fn_apo = sorted(fn_apo, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_apo))

plt.rcParams['figure.figsize'] = (12, 8)
for i, fn in enumerate(fn_apo):
    data_fft = np.loadtxt(fn, comments="!")
    # plt.plot(data_fft[:, 0], data_fft[:, 1], label=labels[i])
    plt.plot(data_fft[:, 0], data_fft[:, 1] / np.max(data_fft[:, 1]),
            label=labels[i])

plt.xlim(-1000, 500)
# plt.xlim(-500000, 500000)
plt.legend()



# %%
# dofft plotten
plt.rcParams['figure.figsize'] = (12, 8)
for i, fn in enumerate(fn_apo):
    data_fft = np.loadtxt(fn, comments="!")
    plt.plot(data_fft[:, 1], data_fft[:, 2], label=labels[i])

plt.xlim(-5000, 5000)
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

    spekfwhm_t1 = spekfwhm - 1/spekt1/np.pi
    return spektemps, spekfwhm, spekt1, spekfwhm_t1


def get_fwhm_temp_t1_2():
    homedir = "/home/karajan/uni/master/analyse/data/crn/data/"
    fwhm_data = np.loadtxt(homedir + "SPEK/spek_342K_380K.data")
    t1_data = np.loadtxt(homedir + "T1/t1_342K_380K.data")
    temp_f = fwhm_data[:,1]
    gamma = fwhm_data[:,5]
    t1 = t1_data[:,3]
    fwhm_data2 = np.loadtxt(homedir + "SPEK/spek_360K_440K.data")
    t1_data2 = np.loadtxt(homedir + "T1/t1_360K_440K_170807.data")
    temp_f2 = fwhm_data2[:,1]
    gamma2 = fwhm_data2[:,5]
    t12 = t1_data2[:,3]

    return np.hstack([temp_f, temp_f2[6:]]), np.hstack([(gamma-1/t1[1::2]/2/np.pi)*2, (gamma2[6:]-1/t12[6:]/2/np.pi)*2])


# # test ob richtige Temperaturen gefunden wurden
# plt.yscale("log")
# plt.scatter(t1temps, t1t1)
# plt.scatter(spektemps, spekt1)



# %%
# T1 rausrechnen
spektemps, spekfwhm, spekt1, spekfwhm_t1 = get_fwhm_temp_t1()
spektemps_2, spekfwhm_t1_2 = get_fwhm_temp_t1_2()
plt.scatter(spektemps, spekfwhm, label="FWHM")
# plt.scatter(spektemps_2, spekfwhm_t1_2, label="FWHM mit $T_1$-Korrektur")

plt.ylim(0, 50000)
plt.xlabel("Temperatur [K]")
plt.ylabel("FWHM [Hz]")
# plt.title("CRN FWHM mit $T_1$-Korrektur")
plt.title("CRN FWHM")
plt.legend(loc=3)

save_plot(plt, "/home/karajan/uni/master/analyse/plots/SPEK/spek_fwhm")


# %%
# Sim FWHM
plt.rcParams['figure.figsize'] = (12, 8)
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
plt.scatter(T_VFT_4er(labels / 128), fwhm*100*2*np.pi, color="r")


# %%
# FWHM plotten
plt.rcParams['figure.figsize'] = (12, 8)


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
        # "/home/karajan/uni/master/analyse/data/bruker/SPEK/bruker_spek.data",
    ]

    for i, file_name in enumerate(file_names):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        gamma = data[:, 5] * 2 / 1e3
        gamma_err = data[:, 6] * 2 / 1e3

        plt.errorbar(temp, gamma, yerr=gamma_err, color="tab:blue", fmt='x')


plot_fwhm()
# plt.xlim(340, 350)
# plt.scatter(T_VFT_4er(labels), fwhms / 125.0, color="g")
# plt.scatter(T_VFT_4er(labels / 128), fwhms / 125.0, color="r")



# %%
# plot Maxima
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

def plot_maxima():
    plt.xlabel("Temperatur [K]")
    # plt.ylabel("FWHM [kHz]")
    plt.ylabel("Maximum [kHz]")
    plt.title("Spektren")

    for i, file_name in enumerate(spek_fn):
        data = np.loadtxt(file_name)
        temp = data[:, 1]
        gamma = data[:, 7] * 2 / 1e3
        gamma_err = data[:, 8] * 2 / 1e3

        plt.errorbar(temp, gamma, yerr=gamma_err, fmt='.')


plot_maxima()

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
            plt.scatter(temps, fwhm / 1e3, marker="o")
        elif kind == "maxim":
            plt.scatter(temps, maxims / 1e3, marker="o")
        elif kind == "mean":
            plt.scatter(temps, means / 1e3, marker="o")


directory = "/home/karajan/uni/master/analyse/data/lineshape/joined"

fn_apo = glob.glob(directory + "/*.fid.dat.spec.nmr")
fn_apo = sorted(fn_apo, key=lambda a: float(a.split("!")[1]))
labels = np.array(sorted(float(fn.split("!")[1]) for fn in fn_apo))

maxims = np.empty(len(fn_apo))
means = np.empty(len(fn_apo))
fwhm = np.empty(len(fn_apo))
for i, fn in enumerate(fn_apo):
    maxims[i], means[i], fwhm[i], temp = get_dofft_info(fn)


print(labels)

plt.rcParams['figure.figsize'] = (12, 8)

factor = 1000
# # FWHM Exp / Sim Vergleich
# plt.scatter(T_VFT_4er(labels), fwhm / factor, color="r")
# # plt.scatter(T_VFT_aug(labels / 128), fwhm/8, color="r")
# # plt.scatter(T_mauro_aug(labels / 128), fwhm/8, color="y")
# plot_fwhm()
# # spektemps, spekfwhm, spekt1, spekfwhm_t1 = get_fwhm_temp_t1()
# # plt.scatter(spektemps, spekfwhm_t1 / 1e3, color="g")
# plt.ylim(0, 50)
# # plot_dofft("fwhm")
# plt.xlabel("Temperatur [K]")
# plt.ylabel("FWHM [kHz]")
# # plt.title("CRN FWHM mit $T_1$-Korrektur")
# plt.title("CRN FWHM Vergleich Simulation")
# # plt.xlim(340, 350)
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/SIM/sim_fwhm")
# plt.show()


# # Maxim Exp / Sim Vergleich
# plt.scatter(T_VFT_4er(labels / 128), maxims/8, color="r", marker="^")
# # plt.scatter(T_VFT_aug(labels / 128), maxims/8, color="g", marker="^")
# # plt.scatter(T_mauro_aug(labels / 128), maxims/8, color="y", marker="^")
# plot_dofft("maxim")

# Mean Exp / Sim Vergleich
plt.scatter(T_VFT_4er(labels), means / factor, color="r", marker="^")
# plt.scatter(T_VFT_aug(labels / 128), means / 15, color="g", marker="^")
# plt.scatter(T_mauro_aug(labels / 128), means / 15, color="y", marker="^")
plot_dofft("mean")
plt.ylim(- 25, 5)
plt.xlabel("Temperatur [K]")
plt.ylabel("Schwerpunkt [kHz]")
# plt.title("CRN FWHM mit $T_1$-Korrektur")
plt.title("CRN Schwerpunkt")
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/SPEK/spek_mean")
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/SIM/sim_mean")



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


home_dir = "/home/karajan/uni/master/analyse/data/"
simfn = "lineshape/joined/lineshape_lifetime!3395656250000.0!.fid.dat.spec.nmr"
sim = np.loadtxt(home_dir + simfn, comments="!")
expfn = "lineshape/exp-sel/dummyfolder/CRN_T2_310K_1504_1.ts.spec.nmr"
exp2 = np.loadtxt(home_dir + expfn, comments="!")

plt.plot(sim[:, 0], sim[:, 1] / np.max(sim[:, 1]) * np.max(exp2[:, 1]))
plt.plot(exp2[:, 0], exp2[:, 1])
plt.xlim(-2e5, 2e5)
save_plot(plt, "/home/karajan/uni/master/analyse/plots/SIM/sim_exp_310")


# %%
# Exp Lineshape Vergleich
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")

plt.rcParams['figure.figsize'] = (12, 12)

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
        plt.plot(data_dofft[:, 0],
                 data_dofft[:, 1] / np.max(data_dofft[:, 1]) - i * 0.2,
                 label=temps[i])

plt.legend(loc='right')
plt.xlim(-100000, 100000)
# plt.show()

plt.xlabel("Frequenz [Hz]")
# plt.ylabel("Schwerpunkt [kHz]")
# plt.title("CRN FWHM mit $T_1$-Korrektur")
plt.title("CRN Lineshape")
# save_plot(plt, "/home/karajan/uni/master/analyse/plots/SPEK/spek_lineshape")
