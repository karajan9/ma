import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq
from scipy import stats
import glob, os


def linear(x, b, m):
    return b + m * x


def t1(x, M0, t1, beta, Moff):
    return M0 * (1 - 2*np.exp(-(x/t1)**beta)) + Moff


def t2(x, M0, t2, beta, Moff):
    return M0 * np.exp(-(2*x/t2)**beta) + Moff


def t2_fest(x, M0, t2):
    return M0 * np.exp(-(2*x/t2)**1.0) + 0.0


def f2(tm, tau, beta, Z, M0, Moff):
    return M0 * ((1 - Z) * np.exp(-(tm/tau)**beta) + Z) * np.exp(-(tm/t1)**beta_t1) + Moff


def gauss(x, a, mu, sigma):
    return a * np.exp(-(x-mu)**2/(2*sigma**2))


def lorentz(x, a, gamma, x0):
    return a/np.pi * gamma / (gamma**2 + (x - x0)**2)


def double_lorentz(x, a1, gamma1, x01, a2, gamma2, x02, xoff):
    return a1/np.pi * gamma1 / (gamma1**2 + (x - x01)**2) + \
           a2/np.pi * gamma2 / (gamma2**2 + (x - x02)**2) + xoff


def phase_finder(phase, cmplx):
    return (cmplx * np.exp(-2 * np.pi * phase/360 * 1j)).imag


def load_spek(fn_nmr=None):
    if fn_nmr is None:
        fn_nmr = glob.glob("*.spec.nmr")[0]
    data_nmr = np.loadtxt(fn_nmr, comments="!")
    freq = data_nmr[:,0]
    real = data_nmr[:,1]
    imag = data_nmr[:,2]
    return freq, real, imag


def load_data():
    fn_nmr = glob.glob("*.dat.nmr")[0]
    data_nmr = np.loadtxt(fn_nmr)
    reptime = data_nmr[:,0]
    real = data_nmr[:,5]
    real_err = data_nmr[:,6]
    imag = data_nmr[:,7]
    imag_err = data_nmr[:,8]
    return reptime, real, real_err, imag, imag_err


def load_FWHMs(directory):
    fn_spek = sorted(glob.glob(directory + "/*/*.spec.nmr"))
    size = len(fn_spek)
    temps = np.empty(size)
    fwhms = np.empty(size)
    index = 0
    for index, fn in enumerate(fn_spek):
        temp, fwhm = load_FWHM(fn)
        temps[index], fwhms[index] = temp, fwhm
    return temps, fwhms


def load_FWHM(fn):
    temp = 0
    fwhm = 0
    with open(fn, "r") as spek_file:
        for line in spek_file:
            if "Sample Temperature" in line:
                temp = float(line.split()[-1])
            elif "FWHM" in line:
                fwhm = float(line.split()[-2])
    return temp, fwhm


def load_spek(fn_spek=None):
    if fn_spek is None:
        if os.path.isdir(fn_spek):
            fn_spek = glob.glob("*.spec.nmr")[0]
    data_spek = np.loadtxt(fn_spek, comments="!")
    freq = data_spek[:,0]
    real = data_spek[:,1]
    imag = data_spek[:,2]
    return freq, real, imag


# the data to be indexed needs to be continuous?
def select_data(start, end, data, indexor=None, remove=False, select_indexor=False):
    mask = []
    if indexor is None:
        mask = (data >= start) & (data < end)
    else:
        mask = (indexor >= start) & (indexor < end)
    if remove is True:
        mask = np.invert(mask)
    if select_indexor is True:
        return data[mask], indexor[mask]
    else:
        return data[mask]


def to_cmplx(real, imag):
    cmplx = np.empty(real.shape, dtype=complex)
    cmplx.real = real
    cmplx.imag = imag
    return cmplx


def get_temperature(directory):
    fn = sorted(glob.glob(directory + "/*.info"))[0]
    with open(fn, "r") as info_file:
        for line in info_file:
            if "Sample Temperature" in line:
                return float(line.split()[-1])
    return 0.0


def get_experiment_number(directory):
    return int(directory.rstrip("/").split("/")[-1].split("_")[0])


def phase_fit(cmplx, start_phase):
    phase = leastsq(phase_finder, start_phase, args=cmplx)[0]
    return phase[0], phase_shift(phase, cmplx)


def phase_shift(phase, cmplx):
    return cmplx * np.exp(-2 * np.pi * phase/360 * 1j)


def sort_values(time, data):
    return np.sort(time), data[np.argsort(time)]


def fit_linear(x, y, sigma=None):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return slope, intercept, std_err


def fit_t1(reptime, cmplx, sigma=None):
    p0 = [80, 1e-4, 0.9, 40]
    p = curve_fit(t1, reptime, cmplx.real, p0=p0, sigma=sigma, absolute_sigma=True)
    p_value, p_cov = p
    p_error = np.sqrt(np.diag(p_cov))
    fit = t1(np.sort(reptime), p_value[0], p_value[1], p_value[2], p_value[3])
    return p_value, p_error, fit


def fit_t2(tau, cmplx, sigma=None):
    p0 = [2500, 0.001, 1.1, 0]
    p = curve_fit(t2, tau, np.abs(cmplx), p0=p0, sigma=sigma, absolute_sigma=True)
    p_value, p_cov = p
    p_error = np.sqrt(np.diag(p_cov))
    fit = t2(np.sort(tau), p_value[0], p_value[1], p_value[2], p_value[3])
    return p_value, p_error, fit


def fit_t2_fest(tau, cmplx, sigma=None):
    p0 = [800, 3e-5]
    p = curve_fit(t2_fest, tau, np.abs(cmplx), p0=p0, sigma=sigma, absolute_sigma=True)
    p_value, p_cov = p
    p_error = np.sqrt(np.diag(p_cov))
    fit = t2(np.sort(tau), p_value[0], p_value[1], p_value[2], p_value[3])
    return p_value, p_error, fit


def fit_f2(tm, cmplx, sigma=None):
    startwerte = [1e-4, 0.6, 0.6, 400, 0]
    p = curve_fit(f2, tm, cmplx.real, p0=startwerte, sigma=sigma, absolute_sigma=True)
    p_value, p_cov = p
    p_error = np.sqrt(np.diag(p_cov))
    tau, beta, Z, M0, Moff = p_value
    fit = f2(tm, tau, beta, Z, M0, Moff)
    return p_value, p_error, fit


def fit_gauss(freq, data, sigma=None):
    startwerte = [900, -4000, 15000]
    p = curve_fit(gauss, freq, data, p0=startwerte, sigma=sigma, absolute_sigma=True)
    p_value, p_cov = p
    p_error = np.sqrt(np.diag(p_cov))
    a, mu, sigma_gauss = p_value
    fit = gauss(freq, a, mu, sigma_gauss)
    return p_value, p_error, fit


def fit_lorentz(freq, data, sigma=None):
    startwerte = [1e8, 10000, -1000]
    p = curve_fit(lorentz, freq, data, p0=startwerte, sigma=sigma, absolute_sigma=True)
    p_value, p_cov = p
    p_error = np.sqrt(np.diag(p_cov))
    a, gamma, x0 = p_value
    fit = lorentz(freq, a, gamma, x0)
    return p_value, p_error, fit


def fit_double_lorentz(freq, data, sigma=None):
    startwerte = [1e8, 10000, -1000, 1e9, 1e4, -1e4, -1000]
    p = curve_fit(double_lorentz, freq, data, p0=startwerte, sigma=sigma, absolute_sigma=True)
    p_value, p_cov = p
    p_error = np.sqrt(np.diag(p_cov))
    a1, gamma1, x01, a2, gamma2, x02, xoff = p_value
    fit = double_lorentz(freq, a1, gamma1, x01, a2, gamma2, x02, xoff)
    return p_value, p_error, fit


def plot_fit(time, fit, p, label=""):
    # plt.text(4e-2, -1000, "$T_1$ = {} ms".format(np.round(p[1] * 1000, decimals=3)), fontsize=10)
    # plt.text(4e-2, -1500, "$\\beta$ = {}".format(np.round(p[2], decimals=3)), fontsize=10)
    plt.plot(time, fit, label=label)


def scatter_data(time, data, label="", marker="x"):
    if np.iscomplexobj(data):
        plt.scatter(time, data.real, color="b", marker=marker, label=label+" real")
        plt.scatter(time, data.imag, color="r", marker=marker, label=label+" imag")
    else:
        plt.scatter(time, data, marker="x", label=label)


def scatter_cmplx(time, cmplx, label=""):
    plt.scatter(time, cmplx.real, color="b", marker="x", label=label)
    plt.scatter(time, cmplx.imag, color="r", marker="x", label=label)


def plot_data(time, data, label="", linestyle="-"):
    plt.plot(time, data, linestyle=linestyle, label=label)


def plot_limits(xmin=None, xmax=None, ymin=None, ymax=None):
    xmin_orig, xmax_orig = plt.xlim()
    ymin_orig, ymax_orig = plt.ylim()
    if xmin is None:
        xmin = xmin_orig
    if xmax is None:
        xmax = xmax_orig
    if ymin is None:
        ymin = ymin_orig
    if ymax is None:
        ymax = ymax_orig
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)


def show_plot():
    plt.legend()
    plt.show()


def save_plot(fn=None):
    plt.legend()
    file_types = ["pdf", "svg"]
    fn_save = "/home/jens/Documents/projekte/crn/170713/plots/diff_cc"
    if fn is not None:
        fn_save = fn
    for file_type in file_types:
        plt.savefig(fn_save + "." + file_type, bbox_inches="tight")
