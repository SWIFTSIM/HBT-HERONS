# Analysing a simulation

To run HBT-HERONS, you always need to provide a path to a parameter file (`<PATH_TO_PARAMETER_FILE>`). Several examples of 
different parameter files are provided in [`examples/configs`](https://github.com/SWIFTSIM/HBT-HERONS/tree/master/examples/configs). A description of all accepted  parameters, their default values, and descriptions can be found in the [runtime parameters](parameter_file.md) page.

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

    Adding the second number is important, since otherwise HBT-HERONS will set `<START_OUTPUT> = MinSnapshotIndex`,
    making it re-start from the first output.