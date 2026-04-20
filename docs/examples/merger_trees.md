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
subhalo_evolution  = catalogue.GetTrackEvolution(TrackId_to_follow)[0]
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

# Load the evolution of all progenitors. Note that loading this in the unsorted catalogue may take a while
disrupted_progenitors_evolution = catalogue.GetTrackEvolution(disrupted_progenitors)

fig, ax1 = plt.subplots(1)
for evolution in disrupted_progenitors_evolution:
    ax1.plot(evolution["Snapshot"], evolution["Mbound"] * catalogue.GetMassUnits_Msunh())
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

# Select orphan subhaloes that sunk.
sunk_progenitors  = catalogue.GetSinkProgenitors(TrackId_to_follow)

# Load the evolution of all progenitors. Note that loading this in the unsorted catalogue may take a while
sunk_progenitors_evolution = catalogue.GetTrackEvolution(sunk_progenitors)
fig, ax1 = plt.subplots(1)
for evolution in sunk_progenitors_evolution:
    ax1.plot(evolution["Snapshot"], evolution["Mbound"] * catalogue.GetMassUnits_Msunh())
ax1.set_xlabel('Output Number')
ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
ax1.set_yscale('log')
plt.show()
```

### All progenitors

In the two examples above, we have seen how to identify ***direct*** subhalo progenitors, classified according to whether they underwent disruption or sinking. A combined list of sunken and disrupted progenitors can be obtained by calling `GetAllProgenitors`. Note that when calling this method, it defaults to returning direct and indirect secondary progenitors of the subhalo of interest. The difference between a direct and indirect secondary progenitor is as follows:

*   **Direct secondary progenitor**: the subhalo merges directly with the main progenitor of the subhalo of interest.
*   **Indirect secondary progenitor**: the subhalo does not merge directly with the main progenitor of the subhalo of interest. Instead, it merges with another subhalo which itself eventually merges with the main progenitor of interest, or with another subhalo that eventually merges with the main branch, etc.

One can always limit the progenitors returned by `GetAllProgenitors` by passing the argument `only_direct_progenitors=True`.

In the code below we show how to identify all direct and indirect secondary progenitors of the subhalo used in the [first example](#main-progenitor).

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

# Select orphan subhaloes that have a connection to the main progenitor brach.
all_progenitors  = catalogue.GetAllProgenitors(TrackId_to_follow)

# Load the evolution of all progenitors. Note that loading this in the unsorted catalogue may take a while
all_progenitor_evolution = catalogue.GetTrackEvolution(all_progenitors)

fig, ax1 = plt.subplots(1)
for evolution in all_progenitor_evolution:
    ax1.plot(evolution["Snapshot"], evolution["Mbound"] * catalogue.GetMassUnits_Msunh())
ax1.set_xlabel('Output Number')
ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
ax1.set_yscale('log')
plt.show()
```

## Merger mass ratios

As shown above, identifying the progenitors of subhaloes and their evolution is relatively simple. This means that statistics that go beyond the evolution of individual subhaloes are correspondingly straightforward to obtain.

In the example below, we show how to measure the cumulative mass ratio distribution of all subhaloes that directly merged with the subhalo used in the [first example](#main-progenitor).

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

# Get the evolution of the main progenitor
evolution_main_progenitor = catalogue.GetTrackEvolution(TrackId_to_follow)[0]

# Select orphan subhaloes that directly merged with our subhalo of interest, and load their peak mass and when it was reached.
all_direct_progenitors = catalogue.GetAllProgenitors(TrackId_to_follow, only_direct_progenitors=True)
all_direct_progenitor_properties = catalogue.LoadSubhalos(subhalo_index = all_direct_progenitors, property_selection=["SnapshotOfLastMaxMass", "LastMaxMass"])

# Compute the ratio between peak mass of each progenitor relative to the bound mass of the main progenitor at the same output.
mass_ratios = - np.ones(len(all_direct_progenitor_properties), float)
for i, progenitor in enumerate(all_direct_progenitor_properties):
    mass_ratios[i] = progenitor["LastMaxMass"] / evolution_main_progenitor["Mbound"][evolution_main_progenitor["Snapshot"] == progenitor["SnapshotOfLastMaxMass"]]

fig, ax1 = plt.subplots(1)
ax1.plot(np.sort(mass_ratios)[::-1], np.arange(len(mass_ratios)) + 1, label="All subhaloes")
ax1.set_ylabel(r"$N(\geq f_{\rm ratio})$")
ax1.set_xlabel(r"$f_{\rm ratio}$")
ax1.set_yscale('log')
ax1.set_xscale('log')
plt.show()
```