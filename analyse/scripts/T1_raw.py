# %%
import os
import sys
home_dir = "/home/karajan/uni/master/ma/analyse"
sys.path.append(os.path.abspath(home_dir + "/scripts"))
from nmr_lib import *

from lmfit import Model


# %%
def analyze_data(directory, verbose=True, normed=False):
    print(directory)
    reptime, real, real_err, imag, imag_err = load_data(directory)
    cmplx = to_cmplx(real, imag)
    phase, cmplx = phase_fit(cmplx, 30)
    reptime, cmplx = sort_values(reptime, cmplx)

    temp = get_temperature(directory)
    experiment_number = get_experiment_number(directory)

    t1_model = Model(t1)
    params = t1_model.make_params()
    params["M0"].set(value=300.0, min=0.0)
    params["t1"].set(value=1e-4, min=1e-7, max=1e-1)
    params["beta"].set(value=1.0, min=0.0, max=2.0)
    params["Moff"].set(value=0.0)
    result = t1_model.fit(cmplx.real, params, x=reptime, weights=1/real_err)

    if data_id == "170807": #  and result.params["beta"].value > 1.99:
        t1_model = Model(t1)
        params = t1_model.make_params()
        params["M0"].set(value=300.0, min=0.0)
        params["t1"].set(value=1e-4, min=1e-7, max=1e-1)
        params["beta"].set(value=1.0, min=0.0, max=2.0, vary=False)
        params["Moff"].set(value=0.0)
        result = t1_model.fit(cmplx.real, params, x=reptime, weights=1/real_err)

    fit = result.best_fit
    if verbose:
        if normed:
            vmax = result.params["M0"].value + result.params["Moff"].value
            vmin = - result.params["M0"].value + result.params["Moff"].value
            cmplx.real = (cmplx.real - vmin) / (vmax-vmin) * 2 - 1
            fit = (result.best_fit - vmin) / (vmax-vmin) * 2 - 1

        plt.scatter(reptime, cmplx.real, label="T = {}".format(np.round(temp, 1)))
        # plt.plot(reptime, cmplx.imag, 'ro')
        plt.plot(reptime, fit)
        plt.xscale("log")
        plt.show()
        print(experiment_number, temp)
        print(result.fit_report())

    return result, experiment_number, temp, phase


def get_analyse(data_dir, verbose=True, resultfile="", normed=False):
    dirs = sorted(glob.glob(data_dir + "/*/"))

    if resultfile != "":
        # os.chdir(data_dir)
        f = open(resultfile, 'w')
        f.write("# Messungs_ID Temperature[K] Phase[degree] T1[s] T1_err Beta "
                + "Beta_err M0 M0_err Moff Moff_err\n")
        print()

    for i, directory in enumerate(dirs):
        result, experiment_number, temp, phase = analyze_data(directory,
                                                              verbose,
                                                              normed)
        if resultfile != "":
            f.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(
                    experiment_number, temp, phase,
                    result.params["t1"].value, result.params["t1"].stderr,
                    result.params["beta"].value, result.params["beta"].stderr,
                    result.params["M0"].value, result.params["M0"].stderr,
                    result.params["Moff"].value, result.params["Moff"].stderr))

    if resultfile != "":
        f.close()


# %%
# data_id = "170912"
data_id = "170912"
data_dir = home_dir + "/data/" + data_id + "/T1"
# data_dir = home_dir + "/data/170706/T1F2"
os.chdir(data_dir)


# %%
home_dir = "/home/karajan/uni/master/ma/analyse"
data_id = "170807"
data_dir = home_dir + "/data/" + data_id + "/T1"
os.chdir(data_dir)

get_analyse(data_dir, verbose=True, normed=False)
# get_analyse(data_dir, verbose=False, resultfile="T1_" + data_id + "_betafest.data")


# %%
data_id = "170912"
data_dir = home_dir + "/data/" + data_id + "/T1"
os.chdir(data_dir)

# plt.rcParams['figure.figsize'] = (12, 8)
get_analyse(data_dir, verbose=True, normed=True)
plt.xlabel("$t_w$ [s]")
plt.ylabel("Magnetisierung (normiert)")
plt.legend(loc=4)
# plt.title("CRN $T_1$")
# save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/T1/t1_roh2")
