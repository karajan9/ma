# %%
import os
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath(home_dir + "/scripts"))
from nmr_lib import *

from lmfit import Model


# %%
data_dir = home_dir + "/data"
t2files = glob.glob(data_dir + "/*/T2/T2_*.data")
t2files

for file in t2files:
    data = np.loadtxt(file)
    # print(data[:,1])
    plt.errorbar(data[:,1], data[:,3], yerr=data[:,4], fmt='o')
plt.yscale("log")
plt.show()

for file in t2files:
    data = np.loadtxt(file)
    # print(data[:,1])
    plt.errorbar(data[:,1], data[:,5], yerr=data[:,6], fmt='o')
plt.show()
