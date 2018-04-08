# %%
import os
import sys
home_dir = "/home/karajan/uni/master/analyse"
sys.path.append(os.path.abspath(home_dir + "/scripts"))
from nmr_lib import *

from lmfit import Model


# %%
data_dir = home_dir + "/data"
t1files = glob.glob(data_dir + "/*/T1/T1_*.data")
t1files

# %%
for file in t1files:
    data = np.loadtxt(file)
    # print(data[:,1])
    plt.errorbar(data[:,1], 1/data[:,3], yerr=data[:,4], fmt='o')
plt.yscale("log")
plt.show()

for file in t1files:
    data = np.loadtxt(file)
    # print(data[:,1])
    plt.errorbar(data[:,1], data[:,5], yerr=data[:,6], fmt='o')
plt.show()
