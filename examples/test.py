# import pandas as pd
# import matplotlib.pyplot as plt
# from cycler import cycler
# import spk

# #### Example using BAHAMAS f_b-M_200 relation at z=0.125 ####

# z = 0.125
# fb_a = 8.44e-05
# fb_pow = 0.275
# k, sup = sup_model(SO=200, z=z, fb_a=fb_a, fb_pow=fb_pow)
# print(k, sup)

# #### Example using binnedf_b-M_200 from BAHAMAS ####

# fig = plt.figure()
# df = pd.read_csv('BAHAMAS_fb_M200.csv')
# print(df)
# z_grouped = df.groupby('z')

# color = plt.cm.viridis(np.linspace(0, 1, len(z_grouped.groups.keys())))
# custom_cycler = iter(cycler(color=color))
# for z, group in z_grouped:
#     c = next(custom_cycler)['color']
#     print(z)
#     print(group)
#     k, sup = sup_model(SO=200, z=z, M_halo=group.M200, fb=group.f_b, k_max=12)
#     plt.plot(k, sup, label='$z=%.2f$' % z, c=c)

# plt.axhline(1, c='k')
# plt.xlim(1e-1, 12)
# plt.ylim(.75, 1.05)
# plt.legend(loc=3, ncol=3, fontsize=10)
# plt.xscale('log')
# plt.xlabel('$k$  $[h/\\mathrm{Mpc}]$')
# plt.ylabel('$P_\\mathrm{hydro}(k)/P_\\mathrm{DM}(k)$')
# plt.tight_layout()
# plt.show()