# Halo substructure and hierarchy

In this page we showcase how to study the population of subhaloes within a single halo, with particular emphasis to the concept of subhalo hierarchy and subhalo parents.

## Satellite peak mass function

Subhaloes that are in the same halo share the value of `HostHaloId`, allowing us to select satellite systems. Each satellite has an assigned hierarchical depth, which encodes how many connections away a subhalo is from the central. Subhaloes that have <br> `Depth = 0`  are centrals, `Depth = 1` are satellites, `Depth = 2` are satellites-of-satellites, etc

In the provided code, we show how to select satellite subhaloes that are in the same halo as the most massive subhalo at the last output of the simulation. We will use this information to plot the satellite peak mass function, which is further subdivided according to the `Depth` of satellites.

```python
import sys
sys.path.append(f"{HBT_HERONS_PATH}/toolbox")
from HBTReader import HBTReader
import matplotlib.pyplot as plt

use_sorted_catalogues = True
catalogue = HBTReader(f"{SORTED_CATALOGUE_BASE_PATH if use_sorted_catalogues else RAW_CATALOGUE_BASE_PATH}",
                    sorted_catalogues=use_sorted_catalogues)

# Get the TrackId of the most massive subhalo. The reader loads the latest
# available snapshot by default and all subhalo properties.
subhaloes = catalogue.LoadSubhaloProperties()
TrackId_to_follow, HostHaloId = subhaloes["TrackId"][subhaloes["Mbound"].argmax()], subhaloes["HostHaloId"][subhaloes["Mbound"].argmax()]

# Select subhaloes in the halo, and remove the central (depth = 0)
subhaloes_in_halo = subhaloes[(subhaloes["HostHaloId"] == HostHaloId) & (subhaloes["Depth"] != 0)]

# Make the plot
fig, ax1 = plt.subplots(1)

# All satellites together
ax1.plot(np.sort(subhaloes_in_halo["LastMaxMass"])[::-1], np.arange(len(subhaloes_in_halo["LastMaxMass"])) + 1,
         label="All satellite subhaloes")

# Do it according to hierarchical depth
for depth_level in range(1, subhaloes_in_halo["Depth"].max() + 1):
    ax1.plot(np.sort(subhaloes_in_halo["LastMaxMass"][subhaloes_in_halo["Depth"] == depth_level])[::-1],
             np.arange(len(subhaloes_in_halo["LastMaxMass"][subhaloes_in_halo["Depth"] == depth_level])) + 1,
             label=f"Depth = {depth_level:d}  subhaloes")

ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_xlabel(r"$M_{\rm peak} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
ax1.set_ylabel(r"$N(\geq M_{\rm peak})$")
ax1.legend()
plt.show()
```

## Hierarchical connections

Every satellite subhalo has `NestedParentTrackId != -1`, which tells us the `TrackId` of the *immediate* parent subhalo of a satellite subhalo. In other words, the `NestedParentTrackId` tells us which `Depth = D` subhalo is the parent of a subhalo of <br> `Depth = D + 1`. Note that the value of `NestedParentTrackId` reflects whether a previously-associated collection of subhaloes were accreted at the same time by a halo, but the subhalo systems may be far apart due to post-infall orbital scattering.

In the following code, we plot the spatial position all subhaloes in the halo that contains the most massive subhalo at the last output of the simulation. We draw connections between subhaloes following the value of `NestedParentTrackId` to indicate connections between subhaloes.

!!! warning

    We have not taken into account periodic boundary conditions when plotting the positions of subhaloes.

```python
import sys
sys.path.append(f"{HBT_HERONS_PATH}/toolbox")
from HBTReader import HBTReader
import matplotlib.pyplot as plt

use_sorted_catalogues = True
catalogue = HBTReader(f"{SORTED_CATALOGUE_BASE_PATH if use_sorted_catalogues else RAW_CATALOGUE_BASE_PATH}",
                    sorted_catalogues=use_sorted_catalogues)

# Required unit information
mass_units = catalogue.GetMassUnits_Msunh()
length_units = catalogue.GetLengthUnits_Mpch()

# The reader loads the latest available snapshot by default and all subhalo properties.
subhaloes = catalogue.LoadSubhaloProperties()

# Load current and peak bound masses of all subhaloes in the simulation
Mbound = subhaloes["Mbound"] * mass_units
Mpeak  = subhaloes["LastMaxMass"] * mass_units

# Load the positions of subhaloes
Centres = subhaloes["ComovingMostBoundPosition"] * length_units

# Information relating to subhalo hierarchy
Depth  = subhaloes["Depth"]
ParentTrackId = subhaloes["NestedParentTrackId"]

# Select the most massive subhalo in the box and its host halo
TargetTrackId = Mbound.argmax()
HostHaloId = subhaloes["HostHaloId"][TargetTrackId]

# Select all subhaloes that share the same host halo (including central)
subhaloes_in_halo = subhaloes["HostHaloId"] == HostHaloId

# Information needed about the subhaloes in this halo
TrackId_subhaloes_in_halo = subhaloes["TrackId"][subhaloes_in_halo]
Mpeak_subhaloes_in_halo = Mpeak[subhaloes_in_halo]
Depth_subhaloes_in_halo = Depth[subhaloes_in_halo]
ParentTrackId_subhaloes_in_halo = ParentTrackId[subhaloes_in_halo]
CentreSatellites = Centres[subhaloes_in_halo]

# Make the plot
fig, ax1 = plt.subplots(1)

ax1.scatter(CentreSatellites[:,0],CentreSatellites[:,1],
            lw=0.3, s=0.01 * Mpeak_subhaloes_in_halo**(1/3), edgecolor='k',c=Depth_subhaloes_in_halo / Depth_subhaloes_in_halo.max(),alpha=0.7)

# Now connect the subhaloes according to hierarchy
for (current_track_id, parent) in zip(TrackId_subhaloes_in_halo,ParentTrackId_subhaloes_in_halo):
    if parent == -1:
        continue
    ax1.plot([CentreSatellites[TrackId_subhaloes_in_halo == current_track_id][0][0],
              CentreSatellites[TrackId_subhaloes_in_halo == parent][0][0]],
             [CentreSatellites[TrackId_subhaloes_in_halo == current_track_id][0][1],
              CentreSatellites[TrackId_subhaloes_in_halo == parent][0][1]],
                lw=0.05,color='k')

ax1.set_xlabel(r"$x \; [\mathrm{Mpc}h^{-1}]$")
ax1.set_ylabel(r"$y \; [\mathrm{Mpc}h^{-1}]$")

plt.show()
```