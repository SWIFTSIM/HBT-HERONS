#!/bin/env python

import matplotlib
matplotlib.use('Agg')

import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt

snap_nr=int(sys.argv[1])
branch=sys.argv[2]

filename = f"/cosma8/data/dp004/jch/failed_mergers/COLIBRE/L0050N0376/Thermal_non_equilibrium/{branch}/failed_mergers/failed_mergers_{snap_nr:03d}.hdf5"
with h5py.File(filename, "r") as f:
    fm = f["failed_merger"][...]
    fm_idx = f["failed_merger_idx"][...]
    r_sep = f["r_sep"][...]
    v_sep = f["v_sep"][...]
    r_sep_parent = f["r_sep_parent"][...]
    v_sep_parent = f["v_sep_parent"][...]

# Also read TrackId etc
merger_file = filename.replace("/failed_mergers_", "/merging_info_")
with h5py.File(merger_file, "r") as f:
    trackid = f["TrackId"][...]
    nbound = f["Nbound"][...]

#
# Plot things which should have merged but didn't
#
plt.figure(figsize=(8,10))
plt.subplot(2,1,1)
plt.plot(r_sep[fm==1], v_sep[fm==1], "r.", label="Not linked by HBT")
plt.plot(r_sep[fm==2], v_sep[fm==2], "g.", label="Linked by HBT")
plt.xscale("log")
plt.xlabel("Normalised position offset")
plt.yscale("log")
plt.ylabel("Normalised velocity offset")
plt.legend(loc="upper right")

# Plot line r_sep + v_sep = 2
r_points = np.logspace(-2.0, 1.0, 100)
v_points = 2 - r_points
plt.plot(r_points, v_points, "k:")

plt.xlim(0.01, 100.0)
plt.ylim(0.01, 100.0)
plt.title("Subhalos which should have merged but didn't")

# Identify subhalos which should have been merged which
# were linked by HBT
weird_halo = (r_sep+v_sep < 2.0) & (fm==2)
i1 = np.arange(len(fm), dtype=int)[weird_halo]
i2 = fm_idx[weird_halo]
for i, j in zip(i1, i2):
    print(i, j, "TrackId = ", trackid[i], trackid[j], "Nbounds = ", nbound[i], nbound[j])

#
# Plot each subhalo's phase space distance to parent
#
plt.subplot(2,1,2)

have_parent = (r_sep_parent >= 0)

plt.plot(r_sep_parent[have_parent], v_sep_parent[have_parent], "k,", rasterized=True)
plt.xscale("log")
plt.xlabel("Normalised position offset")
plt.yscale("log")
plt.ylabel("Normalised velocity offset")

# Plot line r_sep + v_sep = 2
r_points = np.logspace(-2.0, 1.0, 100)
v_points = 2 - r_points
plt.plot(r_points, v_points, "k:")

plt.xlim(0.01, 100.0)
plt.ylim(0.01, 100.0)
plt.title("Phase space distance to parent for all subhalos")

plt.tight_layout()
plt.savefig(f"./plots/failed_mergers_{branch}_{snap_nr:03d}.png")

#plt.show()


