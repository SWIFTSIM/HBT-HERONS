The latest version is available in the [official repository](https://github.com/SWIFTSIM/HBT-HERONS). 
To get a local copy of the master branch, run:

```bash 
git clone https://github.com/SWIFTSIM/HBT-HERONS.git && cd HBT-HERONS
```

As HBT-HERONS saves the commit hash and its status when the code is compiled, it is 
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
    the options relevant to your setup, press `c` (configure) again and then `g` (generate).

### Type of simulation

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :-----------------------------------------------  | :----------------------------------------------- |
| `HBT_DM_ONLY`                    | `Off`| Enable if the simulation is dark matter only.                   |
| `HBT_UNBIND_WITH_THERMAL_ENERGY` | `Off`| Enable to include thermal energy of gas when calculating its binding energy. If enabled, the dataset needs to be loaded from the particle outputs. |


### Internal precision

These flags are used to determine the internal dtypes used by HBT-HERONS, which always
default to using the maximum precision possible. If you would like to reduce the memory footprint
of the code, you can do so using these flags.
