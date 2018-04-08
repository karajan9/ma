# !/usr/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import scipy

plt.rcParams['axes.grid'] = False
plt.rcParams['figure.figsize'] = (25.0 / 2.54, 35.0 / 2.54)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['figure.max_open_warning'] = 200
plt.rcParams['text.usetex'] = True

td_1 = 16384
td_2 = 512
dw = 0.5e-6
in0 = 4e-6
grpdly = 68 + 8

phase_cos = 92.0 * np.pi / 180.0  # Phase for the phase correction
phase_sin = 2.0 * np.pi / 180.0
apodisation = 500
sinweight = 1.0

fftpoints_f1 = 4096
fftpoints_f2 = 1024

spekwidth_min = -70000.0
spekwidth_max = 50000.0
plotfrom = 0e-6
plotto = 400e-6

dw = dw * 2

l11 = 1.0
l12 = 70.0
l21 = 2.0
l22 = 50.0

dshape = (int(td_2 / 2), 2, int(td_1 / 2), 2)

f = open("ser", "rb")
data = np.frombuffer(f.read(), dtype="<i4")
data.shape = dshape

tempdata = np.zeros(dshape)
tempdata[0:48, :, :, :] = data[0:48, :, :, :]

cropped = np.zeros(dshape)
cropped[:, :, 0:int(td_1 / 2) - grpdly, :] = tempdata[:, :, grpdly:, :]

times_1 = dw * np.arange(0, int(td_1 / 2))
times_2 = in0 * np.arange(0, int(td_2 / 2))

mesh_t2, mesh_t1 = np.meshgrid(times_1, times_2)

signal_cc = cropped[:, 0, :, 0] + 1j * cropped[:, 0, :, 1]
signal_ss = cropped[:, 1, :, 0] + 1j * cropped[:, 1, :, 1]

levels_cc = np.linspace(
    np.min(np.array([np.min(np.real(signal_cc)),
                     np.min(np.imag(signal_cc))])),
    np.max(np.array([np.max(np.real(signal_cc)),
                     np.max(np.imag(signal_cc))])), 30)
levels_ss = np.linspace(
    np.min(np.array([np.min(np.real(signal_ss)),
                     np.min(np.imag(signal_ss))])),
    np.max(np.array([np.max(np.real(signal_ss)),
                     np.max(np.imag(signal_ss))])), 30)

prt11 = int(plotfrom / dw)
prt12 = int(plotto / dw)
prt21 = int(plotfrom / in0)
prt22 = int(plotto / in0)

# fitgtimesignalline=plt.figure(dpi=1200)
# ax1=fitgtimesignalline.add_subplot(4,1,1)
# ax2=fitgtimesignalline.add_subplot(4,1,2)
# for i in range(0,8):
# ax1.plot(times_1[prt11:prt12], np.real(signal_cc[i,prt11:prt12]), label=r"$t_p=%.0f$us cos real"%(times_2[i]*1e6))
# ax2.plot(times_1[prt11:prt12], np.imag(signal_cc[i,prt11:prt12]), label=r"$t_p=%.0f$us cos imag"%(times_2[i]*1e6))

# ax1.legend(loc='upper right')
# ax2.legend(loc='upper right')

# ax1=fitgtimesignalline.add_subplot(4,1,3)
# ax2=fitgtimesignalline.add_subplot(4,1,4)
# for i in range(0,8):
# ax1.plot(times_1[prt11:prt12], np.real(signal_ss[i,prt11:prt12]), label=r"$t_p=%.0f$us sin real"%(times_2[i]*1e6))
# ax2.plot(times_1[prt11:prt12], np.imag(signal_ss[i,prt11:prt12]), label=r"$t_p=%.0f$us sin imag"%(times_2[i]*1e6))

# ax1.legend(loc='upper right')
# ax2.legend(loc='upper right')

# fitgtimesignalline.savefig('graphics/signalline.pdf')

# figtimesignals3d=plt.figure(dpi=1200)
# ax1=figtimesignals3d.add_subplot(2,2,1, projection='3d')
# ax1.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_cc[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# ax2=figtimesignals3d.add_subplot(2,2,2, projection='3d')
# ax2.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_cc[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# ax3=figtimesignals3d.add_subplot(2,2,3, projection='3d')
# ax3.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_ss[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# ax4=figtimesignals3d.add_subplot(2,2,4, projection='3d')
# ax4.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_ss[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# figtimesignals3d.savefig('graphics/signal3d.png')

# figtimesignalscont=plt.figure(dpi=1200)
# ax1=figtimesignalscont.add_subplot(2,2,1)
# ax1.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_cc[prt21:prt22,prt11:prt12]), levels_cc, cmap=cm.jet)

# ax2=figtimesignalscont.add_subplot(2,2,2)
# ax2.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_cc[prt21:prt22,prt11:prt12]), levels_cc, cmap=cm.jet)

# ax3=figtimesignalscont.add_subplot(2,2,3)
# ax3.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_ss[prt21:prt22,prt11:prt12]), levels_ss, cmap=cm.jet)

# ax4=figtimesignalscont.add_subplot(2,2,4)
# ax4.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_ss[prt21:prt22,prt11:prt12]), levels_ss, cmap=cm.jet)

# figtimesignalscont.savefig('graphics/signalcontour.pdf')

# Baseline Correction:
for i in range(0, int(td_2 / 2)):
    signal_cc[i, :] = signal_cc[i, :] - signal_cc[i, int(td_1 / 4):].mean()
    signal_ss[i, :] = signal_ss[i, :] - signal_ss[i, int(td_1 / 4):].mean()
for i in range(0, int(td_1 / 2)):
    signal_cc[:, i] = signal_cc[:, i] - signal_cc[int(td_2 / 8 * 3):, i].mean()
    signal_ss[:, i] = signal_ss[:, i] - signal_ss[int(td_2 / 8 * 3):, i].mean()

# Half first Datapoints:
signal_cc[:, 0] = signal_cc[:, 0] / 2.0
signal_cc[0, :] = signal_cc[0, :] / 2.0
signal_ss[:, 0] = signal_ss[:, 0] / 2.0
signal_ss[0, :] = signal_ss[0, :] / 2.0

# Suggestion for Phase correction:
# Coscos:
sugested_cc = -np.angle(signal_cc[0:3, 0:5].mean())
sugested_ss = -np.angle(signal_ss[1:3, 0:5].mean())
print("Phase suggestion for cos-cos: %.0f degree (%.2f rad)" %
      (sugested_cc * 180 / np.pi, sugested_cc))
print("Phase suggestion for sin-sin: %.0f degree (%.2f rad)" %
      (sugested_ss * 180 / np.pi, sugested_ss))

# Phase correction:
signal_cc = signal_cc * np.exp(1j * phase_cos)
signal_ss = signal_ss * np.exp(1j * phase_sin)

# Apodisation:
apo_function = np.exp(-1.0 / 2.0 * (mesh_t2 * apodisation * 2 * np.pi)**2 -
                      1.0 / 2.0 * (mesh_t1 * apodisation * 2 * np.pi)**2)
signal_cc = signal_cc * apo_function
signal_ss = signal_ss * apo_function

# fitgtimesignallinecor=plt.figure(dpi=1200)
# ax1=fitgtimesignallinecor.add_subplot(4,1,1)
# ax2=fitgtimesignallinecor.add_subplot(4,1,2)
# for i in range(0,8):
# ax1.plot(times_1[prt11:prt12], np.real(signal_cc[i,prt11:prt12]), label=r"$t_p=%.0f$us cos real"%(times_2[i]*1e6))
# ax2.plot(times_1[prt11:prt12], np.imag(signal_cc[i,prt11:prt12]), label=r"$t_p=%.0f$us cos imag"%(times_2[i]*1e6))

# ax1.legend(loc='upper right')
# ax2.legend(loc='upper right')

# ax1=fitgtimesignallinecor.add_subplot(4,1,3)
# ax2=fitgtimesignallinecor.add_subplot(4,1,4)
# for i in range(0,8):
# ax1.plot(times_1[prt11:prt12], np.real(signal_ss[i,prt11:prt12]), label=r"$t_p=%.0f$us sin real"%(times_2[i]*1e6))
# ax2.plot(times_1[prt11:prt12], np.imag(signal_ss[i,prt11:prt12]), label=r"$t_p=%.0f$us sin imag"%(times_2[i]*1e6))

# ax1.legend(loc='upper right')
# ax2.legend(loc='upper right')

# fitgtimesignallinecor.savefig('graphics/signallinecor.pdf')

levels_cc = np.linspace(
    np.min(np.array([np.min(np.real(signal_cc)),
                     np.min(np.imag(signal_cc))])),
    np.max(np.array([np.max(np.real(signal_cc)),
                     np.max(np.imag(signal_cc))])), 30)
levels_ss = np.linspace(
    np.min(np.array([np.min(np.real(signal_ss)),
                     np.min(np.imag(signal_ss))])),
    np.max(np.array([np.max(np.real(signal_ss)),
                     np.max(np.imag(signal_ss))])), 30)

# figtimesignals3dcorr=plt.figure(dpi=1200)
# ax1=figtimesignals3dcorr.add_subplot(2,2,1, projection='3d')
# ax1.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_cc[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# ax2=figtimesignals3dcorr.add_subplot(2,2,2, projection='3d')
# ax2.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_cc[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# ax3=figtimesignals3dcorr.add_subplot(2,2,3, projection='3d')
# ax3.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_ss[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# ax4=figtimesignals3dcorr.add_subplot(2,2,4, projection='3d')
# ax4.plot_surface(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_ss[prt21:prt22,prt11:prt12]), ccount=prt12-prt11, rcount=prt22-prt21 , cmap=cm.jet)

# figtimesignals3dcorr.savefig('graphics/signalcorr3d.png')

# figtimesignalscontcorr=plt.figure(dpi=1200)
# ax1=figtimesignalscontcorr.add_subplot(2,2,1)
# ax1.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_cc[prt21:prt22,prt11:prt12]), levels_cc, cmap=cm.jet)

# ax2=figtimesignalscontcorr.add_subplot(2,2,2)
# ax2.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_cc[prt21:prt22,prt11:prt12]), levels_cc, cmap=cm.jet)

# ax3=figtimesignalscontcorr.add_subplot(2,2,3)
# ax3.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.real(signal_ss[prt21:prt22,prt11:prt12]), levels_ss, cmap=cm.jet)

# ax4=figtimesignalscontcorr.add_subplot(2,2,4)
# ax4.contour(mesh_t1[prt21:prt22, prt11:prt12], mesh_t2[prt21:prt22, prt11:prt12], np.imag(signal_ss[prt21:prt22,prt11:prt12]), levels_ss, cmap=cm.jet)

# figtimesignalscontcorr.savefig('graphics/signalcorrcontour.pdf')

# FFT:
# generate frequency mesh:
freq1 = np.fft.fftshift(np.fft.fftfreq(fftpoints_f1, dw))
freq2 = np.fft.fftshift(np.fft.fftfreq(fftpoints_f2, in0))
mesh_f2, mesh_f1 = np.meshgrid(freq1, freq2)

fft_2_cos = np.real(np.fft.fft(signal_cc, n=fftpoints_f1, axis=1))
fft_2_sin = np.imag(np.fft.fft(signal_ss, n=fftpoints_f1, axis=1))

fft_cos = np.real(np.fft.fft(fft_2_cos, n=fftpoints_f2, axis=0))
fft_sin = np.real(np.fft.fft(1j * fft_2_sin, n=fftpoints_f2, axis=0))

fft_sum = fft_cos - sinweight * fft_sin

fft_cos = np.fft.fftshift(fft_cos)
fft_sin = np.fft.fftshift(fft_sin)
fft_sum = np.fft.fftshift(fft_sum)

slice_mask_d1 = np.logical_and(mesh_f2[0, :] >= spekwidth_min,
                               mesh_f2[0, :] < spekwidth_max)
slice_mask_d2 = np.logical_and(mesh_f1[:, 0] >= spekwidth_min,
                               mesh_f1[:, 0] < spekwidth_max)

mesh_f2 = mesh_f2.compress(
    slice_mask_d1, axis=1).compress(
        slice_mask_d2, axis=0)
mesh_f1 = mesh_f1.compress(
    slice_mask_d1, axis=1).compress(
        slice_mask_d2, axis=0)

fft_cos = fft_cos.compress(
    slice_mask_d1, axis=1).compress(
        slice_mask_d2, axis=0)
fft_sin = fft_sin.compress(
    slice_mask_d1, axis=1).compress(
        slice_mask_d2, axis=0)
fft_sum = fft_sum.compress(
    slice_mask_d1, axis=1).compress(
        slice_mask_d2, axis=0)


def symmetrise(spektrum, both=False):
    l1 = len(spektrum)
    l2 = len(spektrum[0])

    if l1 != l2:
        print("Dimensions are not equal! - Dimensions: %i : %i" % (l1, l2))
        return

    for i in range(0, l1):
        for j in range(0, i):
            spektrum[i, j] = (spektrum[i, j] + spektrum[j, i]) / 2.0
            spektrum[j, i] = spektrum[i, j]
            if both:
                spektrum[i, j] = (
                    spektrum[i, j] + spektrum[l1 - j - 1, l2 - i - 1]) / 2.0
                spektrum[l1 - j - 1, l2 - i - 1] = spektrum[i, j]
    return spektrum


fft_cos = symmetrise(fft_cos, both=False)
fft_sin = symmetrise(fft_sin, both=False)
fft_sum = symmetrise(fft_sum, both=False)

randverteilung_cos_d1 = fft_cos.sum(axis=0)
randverteilung_cos_d2 = fft_cos.sum(axis=1)
randverteilung_sin_d1 = fft_sin.sum(axis=0)
randverteilung_sin_d2 = fft_sin.sum(axis=1)
randverteilung_summe_d1 = fft_sum.sum(axis=0)
randverteilung_summe_d2 = fft_sum.sum(axis=1)

cc_max = np.max(fft_cos)
cc_min = np.min(fft_cos)
cc_span = cc_max - cc_min
ss_max = np.max(fft_sin)
ss_min = np.min(fft_sin)
ss_span = ss_max - ss_min
spek_max = np.max(fft_sum)
spek_min = np.min(fft_sum)
spek_span = spek_max - spek_min

fft_cos = (fft_cos - cc_min) / cc_span
fft_sin = (fft_sin - ss_min) / ss_span
fft_sum = (fft_sum - spek_min) / spek_span

cc_max = np.max(fft_cos)
cc_min = np.min(fft_cos)
cc_span = cc_max - cc_min
ss_max = np.max(fft_sin)
ss_min = np.min(fft_sin)
ss_span = ss_max - ss_min
spek_max = np.max(fft_sum)
spek_min = np.min(fft_sum)
spek_span = spek_max - spek_min

levels = np.linspace(cc_min, cc_max - cc_span * 90.0 / 100.0, 40)

fig = plt.figure(figsize=(14.0 / 2.54, 16.0 / 2.54), dpi=1200)

axis1 = fig.add_subplot(2, 2, 1, projection="3d")
axis1.set_ylim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_xlim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_zlim((cc_min, cc_max))
axis1.set_ylim3d((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_xlim3d((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_zlim3d((cc_min, cc_max))
axis1.set_xlabel(r"$ f_1 $ (kHz)")
axis1.set_ylabel(r"$ f_2 $ (kHz)")
axis1.autoscale_view('tight')
pl1 = axis1.plot_surface(
    mesh_f1 / 1000.0,
    mesh_f2 / 1000.0,
    fft_cos,
    rcount=len(fft_cos),
    ccount=len(fft_cos),
    cmap=cm.jet,
    alpha=1.0,
    antialiased=False,
    vmin=spek_min,
    vmax=spek_max)

axis2 = fig.add_subplot(2, 2, 2, aspect='equal')
pl2 = axis2.contour(
    mesh_f1 / 1000.0,
    mesh_f2 / 1000.0,
    fft_cos,
    levels,
    cmap=cm.jet,
    linewidths=[0.2] * len(levels))
axis2.set_ylim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis2.set_xlim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis2.set_xlabel(r"$ f_1 $ (kHz)")
axis2.set_ylabel(r"$ f_2 $ (kHz)")
axis2.autoscale_view('tight')

levels = np.linspace(ss_min, ss_max - ss_span * 90.0 / 100.0, 40)

axis1 = fig.add_subplot(2, 2, 3, projection="3d")
axis1.set_ylim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_xlim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_zlim((ss_min, ss_max))
axis1.set_ylim3d((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_xlim3d((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis1.set_zlim3d((ss_min, ss_max))
axis1.set_xlabel(r"$ f_1 $ (kHz)")
axis1.set_ylabel(r"$ f_2 $ (kHz)")
axis1.autoscale_view('tight')
pl1 = axis1.plot_surface(
    mesh_f1 / 1000.0,
    mesh_f2 / 1000.0,
    fft_sin,
    rcount=len(fft_sin),
    ccount=len(fft_sin),
    cmap=cm.jet,
    alpha=1.0,
    antialiased=False,
    vmin=spek_min,
    vmax=spek_max)

axis2 = fig.add_subplot(2, 2, 4, aspect='equal')
pl2 = axis2.contour(
    mesh_f1 / 1000.0,
    mesh_f2 / 1000.0,
    fft_sin,
    levels,
    cmap=cm.jet,
    linewidths=[0.2] * len(levels))
axis2.set_ylim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis2.set_xlim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axis2.set_xlabel(r"$ f_1 $ (kHz)")
axis2.set_ylabel(r"$ f_2 $ (kHz)")
axis2.autoscale_view('tight')

fig.tight_layout()
fig.subplots_adjust(bottom=0.15)
cbar_ax = fig.add_axes([0.1, 0.08, 0.8, 0.01])
fig.colorbar(pl1, orientation='horizontal', cax=cbar_ax)

fig.text(0.045, 0.94, "cos", size=15)
fig.text(0.045, 0.50, "sin", size=15)

fig.savefig("graphics/dms2dpsektrumcossin.png")

levels = np.linspace(spek_min + spek_span * l11 / 100.0,
                     spek_max - spek_span * l12 / 100.0, 20)

figpaper = plt.figure(figsize=(6.5 / 2.54, 7.5 / 2.54))

axispapermess = figpaper.add_subplot(1, 1, 1, aspect='equal')
axispapermess.set_ylim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axispapermess.set_xlim((spekwidth_min / 1000.0, spekwidth_max / 1000.0))
axispapermess.set_xticks([-70.0, -30.0, 10.0, 50.0])
axispapermess.set_yticks([-70.0, -30.0, 10.0, 50.0])
axispapermess.set_xticklabels(["-70", "-30", "10", "50"])
axispapermess.set_yticklabels(["-70", "-30", "10", "50"])
axispapermess.set_ylabel(r"$ \nu_2 $ (kHz)")
axispapermess.autoscale_view('tight')
plpapermess = axispapermess.contour(
    mesh_f1 / 1000.0,
    mesh_f2 / 1000.0,
    fft_sum,
    levels,
    cmap=cm.jet,
    linewidths=[0.4] * len(levels),
    vmin=0.0,
    vmax=0.5)

figpaper.tight_layout(w_pad=0.5, h_pad=0.5)
figpaper.subplots_adjust(bottom=0.19)
cbar_ax_paper = figpaper.add_axes([0.18, 0.06, 0.77, 0.01])

norm = matplotlib.colors.Normalize(vmin=0.0, vmax=0.5)
sm = plt.cm.ScalarMappable(norm=norm, cmap=plpapermess.cmap)
sm.set_array([])

cb = figpaper.colorbar(sm, orientation='horizontal', cax=cbar_ax_paper)
cb.set_ticks([0.0, 0.25, 0.5])

figpaper.savefig(
    "graphics/dmso2paper1200dpi_apo%.3fl11%.1fl12%.1fl21%.1fl22%.1f.png" %
    (apodisation, l11, l12, l21, l22),
    dpi=1200)

axispapermess.set_ylim((spekwidth_max / 1000.0, spekwidth_min / 1000.0))
axispapermess.set_xlim((spekwidth_max / 1000.0, spekwidth_min / 1000.0))
figpaper.savefig(
    "graphics/dmso2paper1200dpi_inv_apo%.3fl11%.1fl12%.1fl21%.1fl22%.1f.png" %
    (apodisation, l11, l12, l21, l22),
    dpi=1200)
