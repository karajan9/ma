# %%
import os
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath(home_dir + "/scripts"))
from nmr_lib import *

from lmfit import Model


# %%
def analyze_data(directory, t1_data, index, verbose=True):
    print(directory)
    tm, real, real_err, imag, imag_err = load_data(directory)
    cmplx = to_cmplx(real, imag)
    phase, cmplx = phase_fit(cmplx, 30)
    tm, cmplx = sort_values(tm, cmplx)

    temp = get_temperature(directory)
    experiment_number = get_experiment_number(directory)

    if data_id == "170828":
        index = np.argmin(np.abs(t1_data[:, 1] - temp))
        print(index)

    f2t1_model = Model(f2_t1)
    params = f2t1_model.make_params()
    params["M0"].set(value=300.0, min=0.0)
    params["t1"].set(value=t1_data[index, 3], min=1e-7, max=1e-1, vary=True)
    params["beta"].set(value=t1_data[index, 5], min=0.0, max=2.0, vary=True)
    params["Moff"].set(value=0.0)

    # params["M0"].set(value=300.0, min=0.0)
    # params["t1"].set(value=t1_data[3], vary=False)
    # params["beta"].set(value=t1_data[5], vary=False)
    # params["Moff"].set(value=0.0)
    result = f2t1_model.fit(cmplx.real, params, x=tm, weights=1 / real_err)

    if verbose:
        if experiment_number >= 1662:
            pass

        composed = f2_t1(
            tm,
            result.params["M0"].value,
            # np.mean(cmplx.real[0:4]),
            t1_data[index, 3],
            t1_data[index, 5],
            result.params["Moff"].value)
        plt.plot(tm, cmplx.real, 'bo', label="real")
        # plt.plot(tm, cmplx.imag, 'ro', label="imag")
        plt.plot(tm, result.best_fit, 'g-', label="F2 Daten fit")
        plt.plot(
            tm,
            composed,
            'k-',
            label="M0/Moff von fit, $T_1$/$\\beta$ von $T_1$")

        plt.xscale("log")
        plt.xlabel("$t_m$ [s]")
        plt.ylabel("Magnetisierung (a.u.)")
        plt.title("CRN $F_2$")
        plt.legend()

        # save_plot(plt, "/home/karajan/uni/master/analyse/plots/F2/f2_tieftemp")
        plt.show()
        print(experiment_number, temp)
        print(result.fit_report())

    return result, experiment_number, temp, phase, tm, cmplx.real / composed


# %%
def get_analyse(data_dir, verbose=True, resultfile=""):
    dirs = sorted(glob.glob(data_dir + "/*/"))

    if resultfile != "":
        # os.chdir(data_dir)
        f = open(resultfile, 'w')
        f.write("# Messungs_ID Temperature[K] Phase[degree] T1[s] T1_err Beta "
                + "Beta_err M0 M0_err Moff Moff_err\n")

    t1_file = home_dir + "/data/" + data_id + "/T1/T1_" + data_id + ".data"
    if data_id == "170714":
        t1_file = home_dir + "/data/170713/T1/T1_170713.data"
    t1_data = np.loadtxt(t1_file)

    for i, directory in enumerate(dirs):
        if data_id == "170713" or data_id == "170714":
            i *= 2
        result, experiment_number, temp, phase, tm, div = analyze_data(
            directory, t1_data, i, verbose)

        if resultfile != "":
            f.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(
                experiment_number, temp, phase, result.params["t1"].value,
                result.params["t1"].stderr, result.params["beta"].value,
                result.params["beta"].stderr, result.params["M0"].value,
                result.params["M0"].stderr, result.params["Moff"].value,
                result.params["Moff"].stderr))

        div, tm = select_data(0.0, 0.06, div, indexor=tm, select_indexor=True)
        f2t1_model = Model(f2_t1)
        params = f2t1_model.make_params()
        params["M0"].set(value=1.0, min=0.0)
        params["t1"].set(value=1e-3, min=1e-7, max=1e-1, vary=True)
        params["beta"].set(value=1.0, min=0.0, max=2.0, vary=True)
        params["Moff"].set(value=0.0)

        result = f2t1_model.fit(div, params, x=tm)

        plt.xscale("log")
        plt.scatter(tm, div)
        plt.plot(tm, result.best_fit, 'g-', label="F2 Daten fit")
        plt.ylim(-0.1, 1.1)
        plt.show()
        print(result.fit_report())

    if resultfile != "":
        f.close()


# %%
data_id = "170714"
data_dir = home_dir + "/data/" + data_id #+ "/F2SIN"
# data_dir = home_dir + "/data/170706/T1F2"
os.chdir(data_dir)

# %%
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")
get_analyse(data_dir, verbose=True)

# %%
# get_analyse(data_dir, verbose=False, resultfile="T1_" + data_id + "F2.data")

# %%
t1_file = home_dir + "/data/" + data_id + "/T1/T1_" + data_id + ".data"
t1data = np.loadtxt(t1_file)
t1data[0, :]
