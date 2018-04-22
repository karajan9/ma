# %%
import os
import sys
home_dir = "/home/karajan/uni/master/ma/analyse"
sys.path.append(os.path.abspath(home_dir + "/scripts"))
from nmr_lib import *

from lmfit import Model


# %%
def analyze_data(directory, verbose=True, normed=False):
    reptime, real, real_err, imag, imag_err = load_data(directory)
    cmplx = to_cmplx(real, imag)

    if data_id == "170817" or data_id == "170912" or data_id == "170918":
        real_err = select_data(0.0, 2e-5, real_err, indexor=reptime,
                                    remove=True)
        cmplx, reptime = select_data(0.0, 2e-5, cmplx, indexor=reptime,
                                    remove=True, select_indexor=True)

    phase, cmplx = phase_fit(cmplx, 0)
    if np.sum(cmplx.real) < 0:
        cmplx = phase_shift(180, cmplx)
    reptime, cmplx = sort_values(reptime, cmplx)

    temp = get_temperature(directory)
    experiment_number = get_experiment_number(directory)

    t2_model = Model(t2)
    params = t2_model.make_params()
    params["M0"].set(value=300.0, min=0.0)
    params["t2"].set(value=1e-4, min=1e-7, max=1e-1)
    params["beta"].set(value=1.0, min=0.0, max=2.0, vary=False)
    params["Moff"].set(value=0.0)
    result = t2_model.fit(cmplx.real, params, x=reptime, weights=1/real_err)


    x = np.geomspace(1e-6, 5e-2, 100)
    fit = t2_model.eval(result.params, x=x)
    if verbose:
        if normed:
            vmax = result.params["M0"].value + result.params["Moff"].value
            vmin = result.params["Moff"].value
            cmplx.real = (cmplx.real - vmin) / (vmax-vmin)
            fit = (fit - vmin) / (vmax-vmin)

        plt.scatter(
            reptime, cmplx.real, label="T = {}".format(np.round(temp, 1)))
        # plt.plot(reptime, cmplx.imag, 'ro')

        plt.plot(x, fit)
        # plt.plot(reptime, fit)

        plt.xscale("log")
        # plt.show()
        print(experiment_number, temp, phase)
        print(result.fit_report())

    return result, experiment_number, temp, phase


def get_analyse(data_dir, verbose=True, resultfile="", normed=False):
    dirs = sorted(glob.glob(data_dir + "/*/"))

    if resultfile != "":
        # os.chdir(data_dir)
        f = open(resultfile, 'w')
        f.write("# Messungs_ID Temperature[K] Phase[degree] T2[s] T2_err Beta "
                + "Beta_err M0 M0_err Moff Moff_err\n")
        print()

    for i, directory in enumerate(dirs):
        result, experiment_number, temp, phase = analyze_data(directory,
                                                              verbose,
                                                              normed)
        if resultfile != "":
            f.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(
                    experiment_number, temp, phase,
                    result.params["t2"].value, result.params["t2"].stderr,
                    result.params["beta"].value, result.params["beta"].stderr,
                    result.params["M0"].value, result.params["M0"].stderr,
                    result.params["Moff"].value, result.params["Moff"].stderr))

    if resultfile != "":
        f.close()


# %%
home_dir = "/home/karajan/uni/master/ma/analyse"
data_id = "170807"
data_dir = home_dir + "/data/" + data_id + "/T2"
os.chdir(data_dir)

# get_analyse(data_dir, verbose=True, normed=False)
get_analyse(data_dir, verbose=False, resultfile="T2err_" + data_id + "_betafest.data")


# %%
print("."); print("."); print("."); print("."); print("."); print(".")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("."); print("."); print("."); print("."); print("."); print(".")

data_id = "170912"
data_dir = home_dir + "/data/" + data_id + "/T2"
os.chdir(data_dir)

plt.tick_params(which="both", direction="in", top=True, right=True)
get_analyse(data_dir, verbose=True, normed=True)
plt.xlabel("$t_p$ [s]")
plt.ylabel("Magnetisierung (normiert)")
plt.legend(loc=1)
# plt.title("CRN $T_2$")
save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/T2/t2_roh3")

# %%
get_analyse(data_dir, verbose=False, resultfile="T2err_" + data_id + ".data")
