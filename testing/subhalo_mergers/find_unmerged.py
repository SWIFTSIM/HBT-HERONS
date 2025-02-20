#!/bin/env python

import sys

import numpy as np
import h5py


def periodic_distance(pos1, pos2, boxsize):
    dx = pos2 - pos1
    dx[dx >  0.5*boxsize] -= boxsize
    dx[dx < -0.5*boxsize] += boxsize
    return np.sqrt(np.sum(dx**2, axis=-1))


def distance(pos1, pos2, boxsize):
    dx = pos2 - pos1
    return np.sqrt(np.sum(dx**2, axis=-1))


def find_unmerged(input_filename, output_filename):

    # Read in the input catalogue
    with h5py.File(input_filename, "r") as f:
        pos = f["CoreComovingPosition"][...]
        sigma_r = f["CoreComovingSigmaR"][...]
        vel = f["CorePhysicalVelocity"][...]
        sigma_v = f["CorePhysicalSigmaV"][...]
        nbound = f["Nbound"][...]
        boxsize = f["BoxSize"][()]
        trackid = f["TrackId"][...]
        parent_trackid = f["NestedParentTrackId"][...]

    # Build a kdtree from the halo positions
    from scipy.spatial import KDTree
    tree = KDTree(pos, boxsize=boxsize)

    # This will flag any subhalos which should have merged
    failed_merger = np.zeros(len(nbound), dtype=int)
    failed_merger_idx = -np.ones(len(nbound), dtype=int)
    r_sep = -np.ones(len(nbound), dtype=float)
    v_sep = -np.ones(len(nbound), dtype=float)
    r_sep_parent = -np.ones(len(nbound), dtype=float)
    v_sep_parent = -np.ones(len(nbound), dtype=float)

    # Identify the parent halo for each subhalo, where there is one
    import virgo.util.match as m
    parent_index = m.match(parent_trackid, trackid) # will be -1 for top level halos

    min_nbound = 20
    radius_factor = 2.0

    # Loop over subhalos
    for i in range(len(nbound)):

        if nbound[i] < min_nbound:
            continue
        
        # Locate all other subhalos within 2*sigma_r
        max_radius = radius_factor*sigma_r[i]
        ngb_idx = np.asarray(tree.query_ball_point(pos[i,:], max_radius), dtype=int)

        # Remove self and any unresolved neighbours from the list
        keep = (ngb_idx != i) & (nbound[ngb_idx] >= min_nbound)
        ngb_idx = ngb_idx[keep]
    
        if len(ngb_idx) > 0:
            
            # Get normalized phase space distance
            pos_dist = periodic_distance(pos[i,:], pos[ngb_idx,:], boxsize) / sigma_r[i]
            vel_dist = distance(vel[i,:], vel[ngb_idx,:], boxsize) / sigma_v[i]

            # Identify any with velocity difference less than 2*sigma_v
            merged = pos_dist+vel_dist < radius_factor

            # Tag failed mergers: if we satsify the merging condition, the neighbour should have merged
            failed_merger[ngb_idx[merged]] = 1
            failed_merger_idx[ngb_idx[merged]] = i

            # Record separation for problem objects
            r_sep[ngb_idx[merged]] = pos_dist[merged]
            v_sep[ngb_idx[merged]] = vel_dist[merged]

            # Loop over neighbours j which should have merged onto subhalo i but didn't
            for j in ngb_idx[merged]:

                print(f"Subhalo with TrackId={trackid[j]} should have merged to TrackId={trackid[i]} but didn't")
                
                # Check if any of subhalo j's parent subhalos are subhalo i
                k = parent_index[j]
                while k >= 0:
                    print(f"  Parent has TrackId={trackid[k]}")
                    if k == i:
                        failed_merger[j] = 2 # Failed merger despite being in the subhalo hierarchy
                    k = parent_index[k]

        # Get normalized phase space distance to immediate parent halo, if any
        if parent_index[i] >= 0 and nbound[parent_index[i]] > min_nbound:
            r_sep_parent[i] = periodic_distance(pos[i,:], pos[parent_index[i],:], boxsize) / sigma_r[parent_index[i]]
            v_sep_parent[i] = distance(vel[i,:], vel[parent_index[i],:], boxsize) / sigma_v[parent_index[i]]
                    
    # Write out the result
    with h5py.File(output_filename, "w") as outfile:
        outfile["failed_merger"] = failed_merger
        outfile["failed_merger_idx"] = failed_merger_idx
        outfile["r_sep"] = r_sep
        outfile["v_sep"] = v_sep
        outfile["r_sep_parent"] = r_sep_parent
        outfile["v_sep_parent"] = v_sep_parent

if __name__ == "__main__":

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
    find_unmerged(input_filename, output_filename)
