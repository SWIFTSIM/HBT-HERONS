# Analysing a simulation

To run HBT-HERONS, you always need to provide a path to a parameter file (`<PATH_TO_PARAMETER_FILE>`). Several examples of
different parameter files are provided in [`examples/configs`](https://github.com/SWIFTSIM/HBT-HERONS/tree/master/examples/configs). A description of all accepted  parameters, their default values, and descriptions can be found in the [runtime parameters](parameter_file.md) page.

!!! warning

    Given its history-based approach, __HBT-HERONS requires (many) more than one snapshot to produce sensible
    results__. The spacing between outputs should be sufficiently small to
    capture the formation of subhaloes in isolation, which allows the code to tag their particles for subsequent
    identification if they become satellites. If the time spacing is too large, any subhalo that forms and becomes a satellite between two consecutive snapshot outputs will not be extracted
    from the background of the central subhalo that hosts it.

    In the limit of analysing a simulation with a single output, HBT-HERONS becomes a glorified centre finder for Friends-of-Friends groups, as no satellites exist.

## Analyse all outputs

If you want to analyse all of the outputs of the simulation, the following
command will make HBT-HERONS run from the first output (`MinSnapshotIndex`) to the
last one (`MaxSnapshotIndex`) of the simulation:

=== "With MPI"


    ```bash
    mpirun -np <NUMBER_MPI_RANKS> ./HBT <PATH_TO_PARAMETER_FILE>
    ```

=== "Without MPI"

    ```bash
    ./HBT <PATH_TO_PARAMETER_FILE>
    ```
----

You can specify the number of OMP threads to use per MPI rank using `export OMP_NUM_THREADS=<NUMBER_OMP_THREADS>` before
running the command above.

## Analyse a subset of outputs

If you are interested in running HBT-HERONS for a subset of the outputs,
you can specify the range as follows:

=== "With MPI"


    ```bash
    mpirun -np <NUMBER_MPI_RANKS> ./HBT <PATH_TO_PARAMETER_FILE> <START_OUTPUT> <END_OUTPUT>
    ```

=== "Without MPI"

    ```bash
    ./HBT <PATH_TO_PARAMETER_FILE> <PATH_TO_PARAMETER_FILE> <START_OUTPUT> <END_OUTPUT>
    ```
----

where `<START_OUTPUT>` is greater or equal than `MinSnapshotIndex` and `<END_OUTPUT>` is
less or equal to `MaxSnapshotIndex`. This command is also used to continue a run that already
started. In such a case, `<START_OUTPUT>` is the output number
following the last one that was analysed, and `<END_OUTPUT>` equals `MaxSnapshotIndex`.

!!! warning

    Adding a value to `<END_OUTPUT>` is important when restarting, since otherwise HBT-HERONS will set `<END_OUTPUT> = <START_OUTPUT>`, and hence it will only analyse a single output.

!!! warning

    If you re-run HBT-HERONS on a snapshot that is at higher redshift than your latest available HBT-HERONS catalogue, **you will likely need to discard all catalogues that were generated from subsequent snapshots**. *Please verify that you specify the correct snapshots to analyse when continuing an existing HBT-HERONS run.*

    For example, if you have catalogues for snapshots `0 1 2 3 4 5` but restart from snapshot `2`, you might need to remove catalogues `3 4 5`. This happens because the `TrackId` assigned to each subhalo when it is first resolved (`SnapshotIndexOfBirth`) varies from run to run, meaning that subhaloes with `SnapshotIndexOfBirth = 2` will have a different `TrackId` value than from the original run. Thus, the `TrackId` values in catalogues `3 4 5` for those subhaloes will differ from catalogues `2`, making the merger trees based on `TrackId` invalid.