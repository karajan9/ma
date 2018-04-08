import numpy as np
import matplotlib.pyplot as plt
import leastsquarefit

dirnamescos = [
    "35", "34", "38", "37", "36", "53", "39", "55", "49", "51", "45", "57",
    "47"
]
tpcos = [
    1e-6, 5e-6, 10e-6, 15e-6, 20e-6, 25e-6, 30e-6, 35e-6, 40e-6, 50e-6, 60e-6,
    70e-6, 80e-6
]
phasecos = np.array([35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35
                     ]) * np.pi / 180
holdtccos = [
    True, True, True, False, False, False, False, False, False, False, False,
    False, False
]
dirnamessin = [
    "31", "33", "42", "44", "41", "54", "40", "56", "50", "52", "46", "58",
    "48"
]
tpsin = [
    1e-6, 5e-6, 10e-6, 15e-6, 20e-6, 25e-6, 30e-6, 35e-6, 40e-6, 50e-6, 60e-6,
    70e-6, 80e-6
]
phasesin = np.array([35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35
                     ]) * np.pi / 180

trigger_time_f2 = 0e-6

t1measurementfolder = "43"
t1phase = 35 * np.pi / 180
t1trigger_time = 0e-6

#td_1 und td_2 werden mittlerweile automatisch ausgelesen
dw = 1e-6  #Dwelltime, in Topspin als DW bezeichnet, Zeit zwischen zwei Punkten
grpdly = 68  #wie viel ist vor dem eigentlich zu messenden Signal (meist durch Variable grpdly gegeben, welche bei allen euren Messungen 68 war)

T1 = 2.6
beta1 = 1.0
TC = 9.27e-3
betaC = 1.0

maxtplot = 150e-6

dw = dw * 2  #Bruker verdoppelt immer nocheinmal die dwell time


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
    print(
        "Schaetzung der Phase von Ordner %s: %.0f (%.2f rad) verwendet: %.0f" %
        (folder, phase_sug * 180 / np.pi, phase_sug, phase * 180 / np.pi))
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

    #times_2 = times_2[mask]
    #echos=echos[mask]
    #echo_err=echo_err[mask]

    return times_1, signal, times_2, echos, echo_err


def extract_bruker_t1(folder, trigger_time, phase, T1, beta1):
    times_1, signal, times_2, echos, echo_err = extract_bruker_generic(
        folder, trigger_time, phase)

    maxv = np.max(np.real(echos))
    minv = np.min(np.real(echos))
    startparams = [T1, maxv, minv, beta1]
    bounds = ([0.0, -np.inf, -np.inf, 0.0], [np.inf, np.inf, np.inf, 4.0])
    hold = [False, False, False, False]

    fitparams, fiterrs = leastsquarefit.leastsquarefit(
        times_2, np.real(echos), g, startparams, bounds=bounds, hold=hold)

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

    return fitparams, fiterrs


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
        times_2, np.real(echos), fu, startparams, bounds=bounds,
        hold=hold)  # , sigma=np.real(echo_err))

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

    snrerror = np.real(echo_err).mean()

    maxf5v = np.max(np.real(echos[0:4]))
    minf5v = np.min(np.real(echos[0:4]))

    maxm3v = np.max(np.real(echos[14:16]))
    minm3v = np.min(np.real(echos[14:16]))

    maxe3v = np.max(np.real(echos[18:20]))
    mine3v = np.min(np.real(echos[18:20]))

    newerr = np.array([
        snrerror / fitparams[2], (maxm3v - mine3v) / (minf5v - mine3v),
        (minm3v - maxe3v) / (maxf5v - maxe3v),
        ((maxm3v - mine3v) / (minf5v - mine3v) - (minm3v - maxe3v) /
         (maxf5v - maxe3v)) / 2.0
    ])

    return fitparams, fiterrs, newerr


T1fit, T1fiterr = extract_bruker_t1(t1measurementfolder, t1trigger_time,
                                    t1phase, T1, beta1)

T1 = T1fit[0]
beta1 = T1fit[3]

TCcos = []
TCsin = []

newerrfilec = []
newerrfiles = []

for i in range(0, len(dirnamescos)):
    fitpar, fiterr, newerr = extract_bruker_f2(
        dirnamescos[i],
        trigger_time_f2,
        phasecos[i],
        tpcos[i],
        T1,
        beta1,
        TC,
        betaC,
        holdtc=holdtccos[i],
        functionname="cos")
    TCcos.append([
        tpcos[i], fitpar[0], fiterr[0], fitpar[5], fiterr[5], fitpar[4],
        fiterr[4], fitpar[1], fiterr[1], fitpar[6], fiterr[6], fitpar[2],
        fiterr[2], fitpar[3], fiterr[3]
    ])
    newerrfilec.append(newerr)

for i in range(0, len(dirnamessin)):
    fitpar, fiterr, newerr = extract_bruker_f2(
        dirnamessin[i],
        trigger_time_f2,
        phasesin[i],
        tpsin[i],
        T1,
        beta1,
        TC,
        betaC,
        holdtc=False,
        functionname="sin")
    TCsin.append([
        tpsin[i], fitpar[0], fiterr[0], fitpar[5], fiterr[5], fitpar[4],
        fiterr[4], fitpar[1], fiterr[1], fitpar[6], fiterr[6], fitpar[2],
        fiterr[2], fitpar[3], fiterr[3]
    ])
    newerrfiles.append(newerr)

TCcos = np.array(TCcos)
TCsin = np.array(TCsin)

np.savetxt(
    "auswertung/tc_cos.dat",
    TCcos,
    delimiter="\t",
    header="tp TC TCerr betaC betaCerr Z Zerr T1 T1err beta1 beta1err M0 M1")
np.savetxt(
    "auswertung/tc_sin.dat",
    TCsin,
    delimiter="\t",
    header="tp TC TCerr betaC betaCerr Z Zerr T1 T1err beta1 beta1err M0 M1")

np.savetxt("auswertung/newerrc.dat", np.array(newerrfilec), delimiter="\t")
np.savetxt("auswertung/newerrs.dat", np.array(newerrfiles), delimiter="\t")

fig = plt.figure(figsize=(32 / 2.54, 18 / 2.54))

ax1 = fig.add_subplot(3, 1, 1)
ax1.errorbar(TCcos[:, 0], TCcos[:, 1], yerr=TCcos[:, 2], fmt="rx", label="cos")
ax1.errorbar(TCsin[:, 0], TCsin[:, 1], yerr=TCsin[:, 2], fmt="bo", label="sin")
ax1.legend()
ax1.set_yscale("log")
ax1.set_ylim((1e-3, 5e-2))

ax2 = fig.add_subplot(3, 1, 3)
ax2.errorbar(TCcos[:, 0], TCcos[:, 5], yerr=TCcos[:, 6], fmt="rx")
ax2.errorbar(TCsin[:, 0], TCsin[:, 5], yerr=TCsin[:, 6], fmt="bo")
ax2.set_ylim((0, 1))

ax3 = fig.add_subplot(3, 1, 2)
ax3.errorbar(TCcos[:, 0], TCcos[:, 3], yerr=TCcos[:, 4], fmt="rx")
ax3.errorbar(TCsin[:, 0], TCsin[:, 3], yerr=TCsin[:, 4], fmt="bo")
ax3.set_ylim((0, 2))

fig.savefig("auswertung/tc.pdf")
plt.close(fig)
