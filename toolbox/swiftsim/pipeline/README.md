# SWIFT pipeline

Scripts to automate running HBT-HERONS on a SWIFT simulation in COSMA. They generate the
configuration file, and submit the following SLURM jobs, each depending on the previous
one completing successfully:

1. Particle splitting information, if particle splitting was enabled in the simulation
   (see `../particle_splitting`).
2. HBT-HERONS.
3. Sorting of the HBT-HERONS catalogues into a single file per output, with subhaloes in
   ascending TrackId (see `../../catalogue_cleanup`). These are saved in
   `HBT-HERONS/sorted_catalogues`, without particle IDs or potential energies.

# Usage

The scripts must be run from this directory. The simulation folder should contain
`used_parameters.yml` and an empty `HBT-HERONS` folder. To start halo finding:
```bash
./start_halo_finding.sh <PATH_TO_SIMULATION>
```
This copies the submission scripts into the `HBT-HERONS` folder and processes all
snapshots that exist at the time of submission.

To process snapshots that were created after the last submission (e.g. for an ongoing
simulation), run:
```bash
./continue_halo_finding.sh <PATH_TO_SIMULATION>
```
This also submits sorting jobs for any completed HBT-HERONS catalogues which do not have a
sorted counterpart.

The submission scripts in `submission_scripts` assume that a virtual environment was
created by running `../create_cosma_env.sh` within `toolbox/swiftsim`.
