# %%
import os
import sys
home_dir = "/home/karajan/uni/master/ma/analyse"
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

    composed = f2_t1(
        tm,
        result.params["M0"].value,
        # np.mean(cmplx.real[0:4]),
        t1_data[index, 3],
        t1_data[index, 5],
        result.params["Moff"].value)

    if verbose:
        if experiment_number >= 1662:
            pass


        plt.plot(tm, result.best_fit, color="tab:blue")
        plt.scatter(tm, cmplx.real, color="tab:blue", marker="o", label="Messpunkte mit Fit")
        # plt.plot(tm, cmplx.imag, 'ro', label="imag")
        plt.plot(
            tm,
            composed,
            color="tab:orange",
            label="An $F_2$-Daten angepasste $T_1$-Kurve")

        plt.tick_params(which="both", direction="in", top=True, right=True)
        plt.xscale("log")
        plt.xlabel("$t_m$ [s]")
        plt.ylabel("Magnetisierung (a.u.)")
        # plt.title("CRN $F_2$")
        plt.legend(loc=1)

        # save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/F2/f2_310")
        # save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/F2/f2_310")
        plt.show()
        print(experiment_number, temp)
        print(result.fit_report())

    return result, experiment_number, temp, phase, tm, cmplx.real / composed


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

    colors = ["tab:green", "tab:red"]
    labels = [
        "Sin-Sin",
        "Cos-Cos"
    ]
    for i, directory in enumerate(dirs):
        if data_id == "170713" or data_id == "170714":
            index = i * 2
        index = 1
        result, experiment_number, temp, phase, tm, div = analyze_data(
            directory, t1_data, index, verbose)

        if resultfile != "":
            f.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(
                experiment_number, temp, phase, result.params["t1"].value,
                result.params["t1"].stderr, result.params["beta"].value,
                result.params["beta"].stderr, result.params["M0"].value,
                result.params["M0"].stderr, result.params["Moff"].value,
                result.params["Moff"].stderr))

        div2, tm3 = select_data(0.0, 0.06, div, indexor=tm, select_indexor=True, remove=True)
        div, tm2 = select_data(0.0, 0.06, div, indexor=tm, select_indexor=True)
        f2t1_model = Model(f2_t1)
        params = f2t1_model.make_params()
        params["M0"].set(value=1.0, min=0.0)
        params["t1"].set(value=1e-3, min=1e-7, max=1e-1, vary=True)
        params["beta"].set(value=1.0, min=0.0, max=2.0, vary=True)
        params["Moff"].set(value=0.0)

        result = f2t1_model.fit(div, params, x=tm2)

        plt.tick_params(which="both", direction="in", top=True, right=True)
        plt.xscale("log")
        plt.scatter(tm2, div, color=colors[i], label=labels[i])
        plt.scatter(tm3, div2, color=colors[i], marker="x")
        fit = f2_t1(tm, result.params["M0"].value, result.params["t1"].value,
                    result.params["beta"].value, result.params["Moff"].value)
        plt.plot(tm, fit, color=colors[i])
        plt.ylim(-0.1, 1.1)
        plt.legend(loc=3)
        plt.xlabel("$t_m$ [s]")
        plt.ylabel("$M_{F_2} / M_{F_2}'$")
        # plt.show()
        print(result.fit_report())
    save_plot(plt, "/home/karajan/uni/master/ma/analyse/plots/F2/f2_fit3")

    if resultfile != "":
        f.close()


# %%
data_id = "170714"
data_dir = home_dir + "/data/" + data_id #+ "/F2SIN"
# data_dir = home_dir + "/data/170706/T1F2"
os.chdir(data_dir)
get_analyse(data_dir, verbose=False)

# %%
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")
# print("."); print("."); print("."); print("."); print("."); print(".")
data_id = "170710"
data_dir = home_dir + "/data/" + data_id + "/F2COS" + "/290K"
# data_dir = home_dir + "/data/170706/T1F2"
os.chdir(data_dir)
get_analyse(data_dir, verbose=False)


# %%
# get_analyse(data_dir, verbose=False, resultfile="T1_" + data_id + "F2.data")

# %%
t1_file = home_dir + "/data/" + data_id + "/T1/T1_" + data_id + ".data"
t1data = np.loadtxt(t1_file)
t1data[0, :]
