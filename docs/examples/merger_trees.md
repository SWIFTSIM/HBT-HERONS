# Merger tree traversal

In this page we explain how to follow the evolution of a subhalo of interest ([main progenitor](#main-progenitor)), and how to identify subhaloes that contributed to its build-up ([secondary progenitors](#secondary-progenitors)).

## Main progenitor

The `TrackId` of a subhalo is a unique identifier that persists in time throughout the simulation. Thus, following the evolution of a subhalo from the first output when it is found in the simulation (`SnapshotOfBirth`) until the last output when it is resolved as self-bound (`SnapshotOfDeath`) only requires knowing its `TrackId`. In fact, the subhalo can still be tracked as an orphan subhalo after its "death", but only a subset of properties are computed in that case.

In the code below we follow the bound mass evolution of the most massive subhalo, identified at the last output of the simulation.

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
subhaloes = catalogue.LoadSubhalos()
TrackId_to_follow = subhaloes["TrackId"][subhaloes["Mbound"].argmax()]

# Get its bound mass evolution, which returns by default all properties and
# the associated scale factors and snashot output numbers.
subhalo_evolution  = catalogue.GetTrackEvolution(TrackId_to_follow)
Mbound_evolution   = subhalo_evolution["Mbound"] * catalogue.GetMassUnits_Msunh()

fig, ax1 = plt.subplots(1)
ax1.plot(Mbound_evolution)
ax1.set_xlabel('Output Number')
ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
ax1.set_yscale('log')
plt.show()
```

## Secondary progenitors

The catalogues also provide sufficient infomation to identify which subhaloes merged together, helping to establish links between different evolutionary branches through so-called secondary progenitors. Obtaining the secondary progenitors of a given subhalo is more involved that following its main branch. One complication is that subhaloes disappear from the HBT-HERONS catalogues in two different ways:

* [Subhalo disruption](#disrupted-progenitors) as a consequence of subhaloes no longer being self-bound.
* [Subhalo sinking](#sunk-progenitors) when the subhalo is still self-bound but its core coalesces in phase-space with that of another subhalo.

HBT-HERONS provides different information depending on which of the two processes led to the removal of a subhalo from the simulation. Commonly used merger trees do not provide this two-category classification, as following the entire process of
subhalo sinking is difficult and is inadecuately followed in traditional subhalo finders.

### Disrupted progenitors

Disruption occurs when the subhalo is no longer considered to be self-bound. To identify the descendant of a disrupted subhalo, HBT-HERONS uses the `NumTracersForDescendants` most bound tracer particles when the subhalo was last self-bound. The descendant subhalo (`DescendantTrackId`) is the one that is self-bound and contains the majority of the aforementioned particles.

Here we show how to identify all subhaloes that disrupted and merged with the subhalo used in the [first example](#main-progenitor).

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
subhaloes = catalogue.LoadSubhalos()
TrackId_to_follow = subhaloes["TrackId"][subhaloes["Mbound"].argmax()]

# Select orphan subhaloes that disrupted or underwent unresolved sinking.
disrupted_progenitors  = catalogue.GetDisruptionProgenitors(TrackId_to_follow)

# Plot the evolution of a random subset of progenitors
Mbound_evolution = - np.ones((len(disrupted_progenitors), len(catalogue.SnapshotIdList)))
for i, disrupted_progenitor in enumerate(disrupted_progenitors):
    subhalo_evolution = catalogue.GetTrackEvolution(disrupted_progenitor)
    Mbound_evolution[i,subhalo_evolution["Snapshot"]] = subhalo_evolution["Mbound"] * catalogue.GetMassUnits_Msunh()

fig, ax1 = plt.subplots(1)
ax1.plot(Mbound_evolution.T)
ax1.set_xlabel('Output Number')
ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
ax1.set_yscale('log')
plt.show()
```

### Sunk progenitors

The sinking of a subhalo refers to process of its core becoming indistinguishable in phase-space from the core of another subhalo, whilst both remain self-bound. HBT-HERONS identifies a sinking even when two subhaloes are within a threshold distance in phase-space, flagging the least massive of the two as being sunk and converting it to an orphan subhalo. The `TrackId` of the most massive of the two is assigned as the `SinkTrackId` of the subhalo that sunk.

Here we show how to identify all subhaloes that sunk and merged with the subhalo used in the [first example](#main-progenitor).

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
subhaloes = catalogue.LoadSubhalos()
TrackId_to_follow = subhaloes["TrackId"][subhaloes["Mbound"].argmax()]

# Select orphan subhaloes that disrupted or underwent unresolved sinking.
disrupted_progenitors  = catalogue.GetSinkProgenitors(TrackId_to_follow)

# Plot the evolution of a random subset of progenitors
Mbound_evolution = - np.ones((len(disrupted_progenitors), len(catalogue.SnapshotIdList)))
for i, disrupted_progenitor in enumerate(disrupted_progenitors):
    subhalo_evolution = catalogue.GetTrackEvolution(disrupted_progenitor)
    Mbound_evolution[i,subhalo_evolution["Snapshot"]] = subhalo_evolution["Mbound"] * catalogue.GetMassUnits_Msunh()

fig, ax1 = plt.subplots(1)
ax1.plot(Mbound_evolution.T)
ax1.set_xlabel('Output Number')
ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
ax1.set_yscale('log')
plt.show()
```