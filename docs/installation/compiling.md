The latest version is available in the [official repository](https://github.com/SWIFTSIM/HBT-HERONS). 
To get a local copy of the master branch, run:

```bash 
git clone https://github.com/SWIFTSIM/HBT-HERONS.git && cd HBT-HERONS
```

As HBT-HERONS saves the commit hash and whether the code has been modified (see [version tracking](../running/diagnostics.md/#version)) before compiling, it is 
recommended to compile it from its own separate `build` folder. This prevents `git` 
from identifying changes to the files it tracks.
```bash 
mkdir build && cd build 
```

The relevant libraries used to compile the code are found through CMake. Then it is a simple case of compiling to generate the `HBT` executable:

```bash 
cmake ../ && make -j 4
```

Note that there are several important [compile-time options](#compile-time-flags) that can be specified using CMake. Compiling
without specifying their values, as we have just done, uses their default values.


## Compile-time options

We only show the subset of options that are directly relevant to how a simulation
is analysed.

!!! tip

    Using `ccmake` instead of `cmake` opens an interactive terminal view showing
    the all the relevant compile-time options and allows the user to change their
    value more easily. If you are doing this for the first time, the cache will 
    be empty. You need to press `c` (configure) to display the options. After you enable/disable
    the options relevant to your setup, press `c` again and then `g` (generate).

### Type of simulation

Options relating to whether the simulation is hydrodynamical, and whether to include 
the thermal energy of gas during subhalo unbinding.

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :-----------------------------------------------  | :----------------------------------------------- |
| `HBT_DM_ONLY`                    | `Off`| Enable if the simulation is dark matter only.                   |
| `HBT_UNBIND_WITH_THERMAL_ENERGY` | `Off`| Enable to include thermal energy of gas when calculating its binding energy. If enabled, the dataset needs to be loaded from the particle outputs. |

### Internal precision

These options are used to determine the internal data types used by HBT-HERONS, which always
default to using the maximum precision possible. If you would like to reduce the memory footprint
of the code, you can do so using these flags.

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :-----------------------------------------------  | :----------------------------------------------- |
| `HBT_INT8`                     | `Off`| Represent integers using 8 (`On`) or 4  (`Off`) bytes. |
| `HBT_REAL8`                    | `Off`| Represent floats using 8 (`On`) or 4  (`Off`) bytes. |
| `HBT_REAL4_VEL`                    | `Off`| Make the code use 4 bytes for velocities, superseeding the choice made in `HBT_REAL8` if `On`. If `Off`, velocities are represented using the same precision specified by `HBT_REAL8`. |
| `HBT_REAL4_MASS`                    | `Off`| Make the code use 4 bytes for masses, superseeding the choice made in `HBT_REAL8` if `On`. If `Off`, masses are represented using the same precision specified by `HBT_REAL8`.  |
| `HBT_INPUT_INT8`                    | `Off`| Represent particle IDs using 8 (`On`) or 4 (`Off`) bytes. **DEPRECATED**: the code always uses the value specified by `HBT_INT8`. |

### Miscellaneous

Options relating to parallelism and whether to compute the eigenvalues of the 
inertia tensors of each subhalo. Recommended to leave as is.

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :-----------------------------------------------  | :----------------------------------------------- |
| `HBT_USE_OPENMP`                     | `On`| Whether to enable OMP parallelism within the main code. |
| `HBT_USE_GSL`                     | `OFF`| Whether to include the GSL library, which is required to calculate the inertia tensor eigenvalues of each self-bound subhalo. |

### Debugging and testing

Options specifiying parallelism within the unit tests, as well as sanity checks throughout the main code. 

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :-----------------------------------------------  | :----------------------------------------------- |
| `HBT_TEST_MPI_NPROCS`                     | `4`| How many MPI ranks to use within the unit tests. |
| `HBT_TEST_OMP_NUM_THREADS`                    | `4`| How many OMP threads to use within the unit tests. |
| `HBT_CHECK_TRACER_INDEX`                    | `Off`| Explictly do a sanity check to see if the Subhalo.TracerIndex points to the correct particle, which is checked by comparing the number of particles. Important to enable when testing changes that modify the vector of subhalo particles at some point, since the ordering of particles can change and hence TracerIndex be incorrect. |