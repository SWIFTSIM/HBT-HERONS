# Diagnostics

HBT-HERONS prints information used to keep track of its [version](#version) and
[basic progress](#general-output-information) per analysed output. Fine-grained [timing information](#timing-file) 
is also provided in the form of a separate text file.

## Version

The HBT-HERONS executable always prints the branch and commit hash of the code that was used to compile it. 
This information is useful to identify which version of the code was used 
to analyse the simulation, in case it no longer matches the most up to date HBT-HERONS version.
It also states whether the code has been modified without commiting changes (`dirty`) or not (`clean`).

Below is an example based on the master branch from Sep 17, 2024:
```bash
HBT compiled using git branch: master and commit: c31532a40dad661e843abcecadef03144b0cc057 (clean)
```

!!! tip

    To find the corresponding code version in Github for a given commit hash, head to `https://github.com/SWIFTSIM/HBT-HERONS/commit/<COMMIT_HASH>`.

## General output information

Information is printed during the analysis of each simulation output. If something does not look 
right in the output of a HBT-HERONS run, these messages are the first check to see if something has gone wrong (e.g. an unexpected number of particles or Friends of Friends groups). Here is the template that is printed out:

```bash
<TotalNumberParticles> particles loaded at Snapshot <SnapshotIndex> (<SnapshotNumber>)
<TotalNumberFoFGroups> groups loaded at snapshot <SnapshotIndex> (<SnapshotNumber>)
Unbinding...
saving <TotalNumberSubhalos> subhalos to <SubhaloPath>/<SnapshotIndex>
SnapshotIndex <SnapshotIndex> done. It took <TimeToAnalyse> seconds.
```

Each term in brackets refers to the following quantities:

* `<TotalNumberParticles>`: the total number of particles in the simulation output being analysed.
* `<TotalNumberFoFGroups>`: the total number of Friends of Friends groups in the simulation output being analysed.
* `<TotalNumberSubhalos>`: the total number of subhaloes saved in the HBT-HERONS output for that simulation output. 
* `<SnapshotNumber>`: Absolute numbering of the simulation outputs being analysed.
* `<SnapshotIndex>`: Relative numbering of the simulation output being analysed. It is only different from `SnapshotNumber` if a subset of the outputs are being analysed (specified through `SnapshotIdList`). 
* `<SubhaloPath>`: base location where the HBT-HERONS outputs are saved.
* `<TimeToAnalyse>`: total time to analyse the current output, including I/O.

## Timing file 

For a more fine-grained look at how HBT-HERONS spends time within the analysis of each simulation output, a file is saved in `<SubhaloPath>/timing.log`. Each line corresponds to the timing information for a given `SnapshotIndex`. Usually, most of the time is spent on either loading simulation output data or unbinding subhaloes.

The lines are appended to the end of the file, meaning that restarting an analysis will add this information to a pre-existing timing file.