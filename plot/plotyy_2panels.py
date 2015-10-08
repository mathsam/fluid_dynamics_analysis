indx_s=520
indx_e=600

#bt_U = np.mean(mean_u, -1)
fig = plt.figure()
ax1 = fig.add_subplot(211)
y = np.linspace(-np.pi, np.pi, 1024)
ax1_twin = mlab.plotyy(y[indx_s:indx_e], bt_U[indx_s:indx_e], y[indx_s:indx_e], mean_btv2[indx_s:indx_e], None, r'$U_{bt}$', r"$v_{bt}'^2$", ax1)
ax1.legend(loc='upper left')
ax1_twin.legend(loc='upper right', fancybox=False)
ax1_twin.yaxis.grid(True)

ax2 = fig.add_subplot(212)
bt_PV = np.mean(mean_pv, -1)
bt_dPVdy = np.mean(mean_dpvdy, -1)
ax2_twin=mlab.plotyy(y[indx_s:indx_e], bt_PV[indx_s:indx_e], y[indx_s:indx_e], bt_dPVdy[indx_s:indx_e], 'y', r'$PV_{bt}$', r"$dPV_{bt}/dy$", ax2)
ax2.legend(loc='upper left')
ax2_twin.legend(loc='upper right', fancybox=False)
ax2_twin.yaxis.grid(True)
plt.show()