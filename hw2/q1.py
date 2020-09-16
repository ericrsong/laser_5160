import matplotlib.pyplot as plt
import matplotlib as plot
plot.use('Agg')
import numpy as np 
# import helpfunction as hf

# ---------------------- default figure setting ---------------#
# plot.rc('text', usetex=True)
# plot.rc('font', family='cmunrm')
plot.rcParams['mathtext.fontset'] = 'cm'
# plot.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plot.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
plot.rcParams['xtick.major.size'] = 6
plot.rcParams['xtick.major.width'] = 1
plot.rcParams['xtick.minor.size'] = 3
plot.rcParams['xtick.minor.width'] = 1
plot.rcParams['ytick.major.size'] = 6
plot.rcParams['ytick.major.width'] = 1
plot.rcParams['ytick.minor.size'] = 3
plot.rcParams['ytick.minor.width'] = 1
plot.rcParams['axes.linewidth'] = 1.2
plot.rcParams['axes.edgecolor'] = 'black'
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['axes.facecolor'] = 'white'
plot.rcParams['text.color'] = 'black'
plot.rcParams['axes.labelcolor'] = 'black'
plot.rcParams['xtick.color'] = 'black'
plot.rcParams['ytick.color'] = 'black'
# ---------------------- default figure setting ends---------------#

# or_para = hf.hickle_load('C:\\Users\\ys273\\Desktop\\order_parameter_1.h5')
# cid = hf.hickle_load("E:\\vicsek_model_root\\1_uniform_noise\\0_data\\3_joint_plot\\purecid_dict=2_v2.h5")

# eta = np.array([0.01 * i + 0.02 for i in range(362) ])
# eta2 = np.array([0.01 * i + 0.01 for i in range(86) ])
tau1 = 2
tau2 = 1

t = np.array([0.01 * i for i in range(1200)])

N1 = tau1 * (1 + tau1/(tau2 - tau1) * np.exp(-t/tau1) 
        - tau2/(tau2 - tau1) * np.exp(-t/tau2))
N2 = tau2 * (- np.exp(-t/tau2) + 1)


fig, ax1 = plt.subplots()

# color = 'tab:green'
# ax1.axvline(0.53, 0, 1, label=r'$\eta=0.53$', color = 'b')
# ax1.legend()
plt.xlim(0,12)
ax1.set_xlabel(r'$t/\mu s$')
ax1.set_ylabel(r'$N/R_p N_0$') # , color=color)
ax1.plot(t, N1, linestyle='solid', ms=1, color='tab:green', label=r'$N_1$')
ax1.plot(t, N2, linestyle='solid', ms=1, color='tab:red', label=r'$N_2$')
ax1.legend() # loc='center right') 

fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig("test.png")
# ax1.plot(eta2, sf2, 'o', ms=1, color='tab:green', label=r'$q=2$')
# ax1.plot(eta2, sf3, 'o', ms=1, color='tab:brown', label=r'$q=3$')
# ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# color = 'tab:brown'
# ax2.set_ylabel(r"$\textrm{CID}$", color=color)  # we already handled the x-label with ax1
# ax2.plot(eta[0:86], cid[0:86]/np.log2(10), 'o', ms=1, color=color)
# ax2.tick_params(axis='y', labelcolor=color)


# plt.show()
# fig.savefig("..\\2_pipeline\\joint_plot_cid_vs_sf\\tmp\\cid_sf.png")
# fig.savefig("new_order_parameter.png")

plt.cla()

tau1 = 1
tau2 = 2
# t = np.array([0.01 * i for i in range(1000)])

N1 = tau1 * (1 + tau1/(tau2 - tau1) * np.exp(-t/tau1) 
        - tau2/(tau2 - tau1) * np.exp(-t/tau2))
N2 = tau2 * (- np.exp(-t/tau1) + 1)


# fig, ax1 = plt.subplots()

# color = 'tab:green'
# ax1.axvline(0.53, 0, 1, label=r'$\eta=0.53$', color = 'b')
# ax1.legend()
plt.xlim(0,12)
ax1.set_xlabel(r'$t/\mu s$')
ax1.set_ylabel(r'$N/R_p N_0$') # , color=color)
ax1.plot(t, N1, linestyle='solid', ms=1, color='tab:green', label=r'$N_1$')
ax1.plot(t, N2, linestyle='solid', ms=1, color='tab:red', label=r'$N_2$')
ax1.legend() # loc='center right') 
# ax1.plot(eta2, sf2, 'o', ms=1, color='tab:green', label=r'$q=2$')
# ax1.plot(eta2, sf3, 'o', ms=1, color='tab:brown', label=r'$q=3$')
# ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# color = 'tab:brown'
# ax2.set_ylabel(r"$\textrm{CID}$", color=color)  # we already handled the x-label with ax1
# ax2.plot(eta[0:86], cid[0:86]/np.log2(10), 'o', ms=1, color=color)
# ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout() 
fig.savefig("test2.png")