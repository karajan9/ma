import numpy as np
import matplotlib.pyplot as plt
import os

plt.style.use("ggplot")
# plt.rc("font", **{"family":"sans-serif", "sans-serif":["Roboto"]})

directory = "/home/jens/Documents/projekte/crn/170706/T1"
os.chdir(directory)

# fn = "1238_CRN_T1_300K/CRN_T1_300K_1238_1.ts"  # phase & echo_time 16000
fn = "1279_CRN_T1_270K/CRN_T1_270K_1279_20.ts"

data = np.loadtxt(fn)
time = data[:,0]
real = data[:,1]
imag = data[:,3]

start_time = 0e-6
start_index = np.argmin(np.abs(time - start_time))
end_time = 35e-6
end_index = np.argmin(np.abs(time - end_time))
# echo_time = 30e-6
# phase = 0

cmplx = np.empty(real.shape, dtype=complex)
cmplx.real = real
cmplx.imag = imag
cmplx_sel = cmplx[start_index:end_index]

max_value = 0
max_phase = 0
max_time = 0

for ph in np.linspace(0, 360, 200):
    temp = cmplx_sel * np.exp(-2*np.pi * ph/360 * 1j)
    temp_max = np.max(np.abs(temp.real))
    # print(temp_max, ph)
    if temp_max > max_value:
        max_value = temp_max
        max_phase = ph
        max_time = time[np.argmax(np.abs(temp.real)) + start_index]

phase = max_phase
print("Max Value = %f" % max_value)
print("Max Phase = %f" % max_phase)
print("Max Time = %f" % max_time)


cmplx *= np.exp(-2*np.pi * phase/360 * 1j)

# plt.xscale("log")
plt.xlim((0, 1e-4))
# plt.ylim((-80, 80))

plt.plot(time, real, label="real", color="b", linestyle="dashed")
plt.plot(time, imag, label="imag", color="r", linestyle="dashed")
plt.plot(time, cmplx.real, label="shifted real", color="b")
plt.plot(time, cmplx.imag, label="shifted imag", color="r")
# plt.plot(time[start_index:end_index], cmplx.real[start_index:end_index], label="shifted real")
# plt.plot(time[start_index:end_index], cmplx.imag[start_index:end_index], label="shifted imag")

plt.legend()
plt.title("CRN Phasenbestimmung")

plt.show()
