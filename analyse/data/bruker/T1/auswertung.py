import numpy as np
import matplotlib.pyplot as plt
import leastsquarefit
from nmr_lib import *
import scipy.optimize


dirnamet1 = [
    "12", "14", "17", "20", "23", "25", "28", "30", "33", "35", "38", "40",
    "43", "46", "49", "56"
]
tempt1 = [
    380.0, 390.0, 375.0, 385.0, 370.0, 365.0, 360.0, 355.0, 350.0, 345.0,
    340.0, 330.0, 335.0, 320.0, 310.0, 340.0
]

trigger_time_t1 = 4e-6
t1phase = 12 * np.pi / 180

#td_1 und td_2 werden mittlerweile automatisch ausgelesen
dw = 1e-6  # Dwelltime, in Topspin als DW bezeichnet, Zeit zwischen zwei Punkten
grpdly = 68  # wie viel ist vor dem eigentlich zu messenden Signal
# (meist durch Variable grpdly gegeben, welche bei allen euren
# Messungen 68 war)

maxtplot = 150e-6

dw = dw * 2  # Bruker verdoppelt immer nocheinmal die dwell time


#Fit Funktionen:
def fu(x, *param):
    return param[2] * (
        (1 - param[4]) * np.exp(-(x / param[0])**param[5]) + param[4]
    ) * np.exp(-(x / param[1])**param[6]) + param[3]


def g(x, *param):
    return param[1] * np.exp(-(x / param[0])**param[3]) + param[2]


def extract_bruker_generic(folder, trigger_time, phase):
    effgrpdly = int(grpdly + trigger_time / dw)
    if effgrpdly < 0:
        effgrpdly = 0
    times_2 = np.loadtxt("%s/vdlist" % (folder))
    f = open("%s/ser" % (folder), "rb")
    data = np.frombuffer(f.read(), dtype="<i4")
    td_2 = len(times_2)
    td_1 = int(len(data) / td_2)
    dshape = (td_2, int(td_1 / 2), 2)
    data.shape = dshape
    cropped = np.zeros(dshape)
    cropped[:, 0:int(td_1 / 2) - effgrpdly, :] = data[:, effgrpdly:, :]
    signal = cropped[:, :, 0] + 1j * cropped[:, :, 1]
    times_1 = dw * np.arange(0, int(td_1 / 2))

    phase_sug = np.angle(
        signal[0:5, 0:5].mean()
    )  #die gewünschte Phase kann ganz gut abgeschätzt werden indem der Winkel in der komplexen Ebene für den ersten Wert berechnet wird.
    # print(
    #     "Schaetzung der Phase von Ordner %s: %.0f (%.2f rad) verwendet: %.0f" %
    #     (folder, phase_sug * 180 / np.pi, phase_sug, phase * 180 / np.pi))
    #Phasenkorrektur anwenden:
    signal = signal * np.exp(-1j * phase)

    fig = plt.figure(figsize=(32 / 2.54, 18 / 2.54))
    ax = fig.add_subplot(1, 1, 1)

    for i in range(0, len(signal)):
        ax.plot(times_1[0:int(maxtplot / dw)],
                np.real(signal[i, 0:int(maxtplot / dw)]))
        ax.plot(times_1[0:int(maxtplot / dw)],
                np.imag(signal[i, 0:int(maxtplot / dw)]))

    fig.savefig("auswertung/%s_signal.pdf" % (folder))
    plt.close(fig)

    #echos ausrechnen:
    echos = signal[:, 0:10].mean(axis=1)
    echo_err = np.real(signal[:, -513:-1]).std(axis=1) + 1j * np.real(
        signal[:, -513:-1]).std(axis=1)

    mask = echo_err != 0.0

    times_2 = times_2[mask]
    echos = echos[mask]
    echo_err = echo_err[mask]

    return times_1, signal, times_2, echos, echo_err


def extract_bruker_t1(folder, trigger_time, phase, T1, beta1):
    times_1, signal, times_2, echos, echo_err = extract_bruker_generic(
        folder, trigger_time, phase)

    # Phasenkorrektur
    phase, echos = phase_fit(echos, phase)

    maxv = np.max(np.real(echos))
    minv = np.min(np.real(echos))
    # print(maxv, minv)

    startparams = [T1, minv, maxv, beta1]
    bounds = ([1e-7, -1e7, -1e6, 0.0], [1e-2, 1e6, 1e6, 2.0])
    # bounds = ([1e-7, maxv - np.abs(maxv)*0.5, minv - np.abs(minv)*0.5, 0.0],
    #           [1e-2, maxv + np.abs(maxv)*0.5, minv + np.abs(minv)*0.5, 2.0])
    hold = [False, False, False, False]

    fitparams, fiterrs = leastsquarefit.leastsquarefit(
        times_2, np.real(echos), g, startparams, bounds=bounds, hold=hold)

    # fitparams, fiterrs = scipy.optimize.curve_fit(

    # print(fitparams)

    fitcurvet = np.logspace(
        np.floor(np.log10(np.min(times_2))),
        np.ceil(np.log10(np.max(times_2))),
        300,
        endpoint=True)
    fitcurvem = g(fitcurvet, *fitparams)

    np.savetxt(
        "auswertung/%s_T1.fit.dat" % (folder),
        np.array([
            fitcurvet, fitcurvem, fitcurvem / fitparams[1],
            (fitcurvem - fitparams[2]) / fitparams[1]
        ]).transpose(),
        header="t_m M M/M0 (M-M1)/M0")
    np.savetxt(
        "auswertung/%s_T1.dat" % (folder),
        np.array([
            times_2,
            np.real(echos),
            np.real(echo_err),
            np.imag(echos),
            np.imag(echo_err),
            np.real(echos) / fitparams[2],
            np.real(echo_err) / fitparams[2],
            np.imag(echos) / fitparams[2],
            np.imag(echo_err) / fitparams[2]
        ]).transpose(),
        header=
        "t_m Re(M) Re(M_err) Im(M) Im(Merr) Re(M)/M0 Re(Merr)/M0 Im(M)/M0 Im(Merr)/M0"
    )

    fig = plt.figure(figsize=(32 / 2.54, 18 / 2.54))
    ax = fig.add_subplot(1, 1, 1)
    ax.errorbar(times_2, np.real(echos), yerr=np.real(echo_err), fmt="bo")
    ax.errorbar(times_2, np.imag(echos), yerr=np.imag(echo_err), fmt="rx")
    ax.plot(fitcurvet, fitcurvem, "b-")
    ax.set_xscale("log")

    plt.figtext(0.01, 0.05, "Fit Parameters: %s" % np.array_str(fitparams))
    plt.figtext(0.01, 0.01, "Fit Uncertainties: %s" % np.array_str(fiterrs))

    fig.savefig("auswertung/%s_T1.plot.pdf" % (folder))
    plt.close(fig)

    return fitparams, fiterrs, phase


def extract_bruker_f2(folder,
                      trigger_time,
                      phase,
                      tp,
                      T1,
                      beta1,
                      TC,
                      betaC,
                      holdtc=False,
                      functionname=""):
    times_1, signal, times_2, echos, echo_err = extract_bruker_generic(
        folder, trigger_time, phase)

    maxv = np.max(np.real(echos))
    minv = np.min(np.real(echos))
    startparams = [TC, T1, maxv, minv, 0.5, betaC, beta1]
    bounds = ([0.0, 0.0, -np.inf, -np.inf, 0.0, 0.0, 0.0],
              [np.inf, np.inf, np.inf, np.inf, 1.0, 4.0, 4.0])
    hold = [holdtc, True, False, False, False, True, True]

    fitparams, fiterrs = leastsquarefit.leastsquarefit(
        times_2, np.real(echos), fu, startparams,
        bounds=bounds, hold=hold)  # , sigma=np.real(echo_err))

    fitcurvet = np.logspace(
        np.floor(np.log10(np.min(times_2))),
        np.ceil(np.log10(np.max(times_2))),
        300,
        endpoint=True)
    fitcurvem = fu(fitcurvet, *fitparams)

    np.savetxt(
        "auswertung/%s_tp%i_fit_%s.dat" % (folder, int(tp * 1e6),
                                           functionname),
        np.array(
            [fitcurvet, fitcurvem,
             (fitcurvem - fitparams[3]) / fitparams[2]]).transpose(),
        header="t_m M (M-M1)/M0")
    np.savetxt(
        "auswertung/%s_tp%i_%s.dat" % (folder, int(tp * 1e6), functionname),
        np.array([
            times_2,
            np.real(echos),
            np.real(echo_err),
            np.imag(echos),
            np.imag(echo_err),
            np.real(echos - fitparams[3]) / fitparams[2],
            np.real(echo_err) / fitparams[2],
            np.imag(echos - fitparams[3]) / fitparams[2],
            np.imag(echo_err) / fitparams[2]
        ]).transpose(),
        header=
        "t_m Re(M) Re(M_err) Im(M) Im(Merr) Re(M-M1)/M0 Re(Merr)/M0 Im(M-M1)/M0 Im(Merr)/M0"
    )

    fig = plt.figure(figsize=(32 / 2.54, 18 / 2.54))
    ax = fig.add_subplot(1, 1, 1)
    ax.errorbar(times_2, np.real(echos), yerr=np.real(echo_err), fmt="bo")
    ax.errorbar(times_2, np.imag(echos), yerr=np.imag(echo_err), fmt="rx")
    ax.plot(fitcurvet, fitcurvem, "b-")
    ax.set_xscale("log")

    plt.figtext(0.01, 0.05, "Fit Parameters: %s" % np.array_str(fitparams))
    plt.figtext(0.01, 0.01, "Fit Uncertainties: %s" % np.array_str(fiterrs))

    fig.savefig("auswertung/%s_tp%i.plot.pdf" % (folder, int(tp * 1e6)))
    plt.close(fig)

    return fitparams, fiterrs


T1 = []
T1err = []
beta1 = []
beta1err = []

print("# Messungs_ID Temperature[K] Phase[degree] T1[s] T1_error Beta Beta_error M0 M0_error Moff Moff_err")
for i, fold in enumerate(dirnamet1):
    fitparam, fiterrs, phase = extract_bruker_t1(fold, trigger_time_t1, t1phase,
                                          1.0e-3, 1.0)
    T1.append(fitparam[0])
    T1err.append(fiterrs[0])
    beta1.append(fitparam[3])
    beta1err.append(fiterrs[3])

    print("{} {} {} {} {} {} {} {} {} {} {}".format(fold, tempt1[i], phase, fitparam[0],
                    fiterrs[0], fitparam[3], fiterrs[3], 0.0, 0.0, 0.0, 0.0))
    

T1 = np.array(T1)
T1err = np.array(T1err)
beta1 = np.array(beta1)
beta1err = np.array(beta1err)
tempt1 = np.array(tempt1)

fig = plt.figure(figsize=(32 / 2.54, 18 / 2.54))

ax1 = fig.add_subplot(1, 1, 1)
ax1.errorbar(1000.0 / tempt1, T1, yerr=T1err, fmt="rx", label="T1")
ax1.legend()
ax1.set_yscale("log")
ax1.set_ylim((1e-5, 1e-1))

# print(T1)

fig.savefig("auswertung/t1.pdf")
plt.close(fig)
