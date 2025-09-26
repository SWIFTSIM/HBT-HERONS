# Diagnostics

HBT-HERONS prints information used to keep track of its [version](#version) and [basic progress](#logged-information) per analysed output. More fine-grained [timing information](../outputs/timing.md) is provided in the form of a separate text file, and optionally, as a timestamp for each individual subhalo within the subhalo catalogues.

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

## Logged information

Information is printed during the analysis of each simulation output. If something does not look
right in the output of a HBT-HERONS run, these messages are the first check to see if something has gone wrong (e.g. an unexpected number of particles or Friends of Friends groups). Here is the template that is printed out and a corresponding example:

=== "Output example"

    ```bash
    33223356 particles loaded from snapshot 44 (SnapshotIndex = 5)
    18150 FOF groups loaded from snapshot 44 (SnapshotIndex = 5). Largest FOF group has 6621 particles.
        Number of hostless subhaloes = 84
        Number of FOF groups without pre-existing subhaloes = 10303
    Reassigning particles...
        Number of reassigned particles = 15123
    Unbinding... Took 3.1258 seconds. Maximum imbalance across ranks was 1.2778.
        Number of merged subhaloes = 116
        Number of disrupted subhaloes = 637
        Number of newly identified subhaloes = 10047
        Number of FOF groups without any self-bound subhaloes = 256
    Finding descendants of merged and disrupted subhaloes...
        Number of subhaloes assigned a descendant = 739
        Number of subhaloes without identified descendants = 14
    Saving 20952 subhalos to /cosma7/data/dp004/dc-foro1/colibre/HBT-HERONS/testing/test_output/044
    Snapshot 44 (SnapshotIndex = 5) done. It took 27.2268 seconds.

    ```

=== "Output template"

    ```bash
    <NumberParticles> particles loaded from snapshot <SnapshotNumber> (SnapshotIndex = <SnapshotIndex>)
    <NumberFOFGroups> FOF groups loaded from snapshot <SnapshotNumber> (SnapshotIndex = <SnapshotIndex>). Largest FOF group has <NumberOfParticlesLargestFOF> particles.
        Number of hostless subhaloes = <NumberHostlessSubhaloes>
        Number of FOF groups without pre-existing subhaloes = <NumberNewGroups>
    Reassigning particles...
        Number of reassigned particles = <NumberReattachments>
    Unbinding... Took <TimeToUnbind> seconds. Maximum imbalance across ranks was <Imbalance>.
        Number of merged subhaloes = <NumberMergers>
        Number of disrupted subhaloes = <NumberDisruptions>
        Number of newly identified subhaloes = <NumberNewSubhalos>
        Number of FOF groups without any self-bound subhaloes = <NumberUnresolvedGroups>
    Finding descendants of merged and disrupted subhaloes...
        Number of subhaloes assigned a descendant = <NumberSubhalosWithDescendants>
        Number of subhaloes without identified descendants = <NumberSubhalosWithoutDescendants>
    Saving <NumberSubhalos> subhalos to <SubhaloPath>/<SnapshotNumber>
    Snapshot <SnapshotNumber> (SnapshotIndex = <SnapshotIndex>) done. It took <TimeToAnalyse> seconds.
    ```

Each term is intended to provide an idea of how the run is performing, where the outputs are being saved, and some general statistics about FOF groups, subhaloes and merger trees. The meaning of each term is explained below.

### Paths and progress

`SubhaloPath`
:   Base location where the HBT-HERONS outputs are saved.

`SnapshotNumber`
:   Absolute numbering of the simulation outputs being analysed.

`SnapshotIndex`
:   Relative numbering of the simulation output being analysed. It is only different from `SnapshotNumber` if a subset of the simulation outputs are being analysed (specified through the runtime parameter `SnapshotIdList`).

### Particle counts

`NumberParticles`

:   The number of particles of all types in the current simulation output.

`NumberOfParticlesLargestFOF`

:   The number of particles of all types in the largest FOF group (by particle number) in the current simulation output.

`NumberReattachments`
:   The number of particles that have been removed from the source subgroup of their associated subhalo and moved to the source of a different subhalo. Only prints out if `ReassignParticles  = 1`.

### Number of structures

`NumberFOFGroups`
:   The number of Friends of Friends groups in the current simulation output.

`NumberSubhalos`
:   The total number of subhaloes in the current simulation output. Since subhaloes that are no longer self-bound (orphan subhaloes) are also saved, this number reflects to the total number of subhaloes that have ever existed in the simulation.

### Timing and imbalance

`TimeToAnalyse`
:   The time it took to analyse the current simulation output in seconds, including I/O.

`TimeToUnbind`
:   The time it took to unbind all subhaloes of the current simulation output in seconds.

`Imbalance`
:   The ratio between the time of the MPI rank that took the longest to do the unbinding, relative to the average time it took across all ranks. Larger values are worse.

### FOF host finding

`NumberHostlessSubhaloes`
:   The number of subhaloes which do not have an assigned host FOF group (`HostHaloId = -1`) in the current simulation output.

`NumberNewGroups`
:   The number of FOF groups which do not have any pre-existing (resolved or disrupted) subhalo in the current simulation output. These FOF groups will be the seeds for new subhaloes in the catalogues, if those seeds are confirmed to be self-bound.

### Subhalo statistics

`NumberUnresolvedGroups`
: The number of FOF groups that do not contain any resolved subhaloes. This number includes two types of FOF group populations. One type corresponds to FOF groups that only have orphan subhaloes. The second type corresponds to FOF groups that do not have any orphan subhalo, i.e. they were a potential seed for new subhaloes but their particles were not self-bound.

`NumberMergers`
:   The total number of subhaloes the current simulation output that were manually merged by HBT-HERONS, after identifying a phase-space overlap with another subhalo in its hierarchy.

`NumberDisruptions`
:   The total number of subhaloes the current simulation output whose number of bound total (or tracer) particles was below a chosen threshold.

`NumberNewSubhalos`
:   The total number of newly identified subhaloes in the current simulation output.

### Merger tree

`NumberSubhalosWithDescendants`
:   The number subhaloes that disrupted or merged in the current simulation output which have an assigned descendant, i.e. the most bound tracer particles from the previous output are now bound to a different subhalo in the current output.

`NumberSubhalosWithoutDescendants`
:   The number subhaloes that disrupted or merged in the current simulation output which do not have an assigned descendant, i.e. the most bound tracer particles from the previous output are not bound to any subhalo in the current output.