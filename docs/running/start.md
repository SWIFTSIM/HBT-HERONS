# Start analysis

To run HBT-HERONS, you need to provide a path to a parameter file (`<PATH_TO_PARAMETER_FILE>`). Several examples of 
different parameter files are provided in [`examples/configs`](https://github.com/SWIFTSIM/HBT-HERONS/tree/master/examples/configs). A description of all accepted  parameters, their default values, and descriptions can be found in the [runtime parameters](#runtime-parameters) page.

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
running the above command.

!!! warning 

    The simulation outputs need to be sufficiently finely spaced for HBT-HERONS
    to work well. Around 64 outputs between z = 20 and z = 0 works well based on
    testing using the FLAMINGO simulations.