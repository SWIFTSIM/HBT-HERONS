# Subhalo hierarchy

HBT-HERONS organises subhaloes within haloes into a hierarchical tree, with connections that determine how mass is transferred between subhaloes. There are two types of subhaloes, which differ in how they can grow in mass:

*   **Central subhalo**: there is only one central subhalo in a halo at any given time. It is the most massive subhalo and the only one that can acquire diffuse mass recently accreted by the halo. The other two growth modes include [sinking](./subhalo_sinking.md) and [stripping mass](./unbinding.md) from satellite subhaloes.
*   **Satellite subhaloes**: all of the other subhaloes in the halo. They may be directly connected to the central subhalo, or to other satellite subhaloes if they were associated to them prior to their infall to the current halo. Contrary to central subhaloes, they cannow grow by diffuse mass accretion, but they can still become more massive through sinking and stripping mass from other satellite subhaloes in their hierarchical tree branch.

The hierarchy in HBT-HERONS is built using the past evolutionary history. When a new subhalo appears in the simulation, it is always a central subhalo. However, neighbouring haloes may coalesce into a single halo with time, meaning that central subhaloes become satellite subhaloes. Identifying which subhaloes are centrals and which are satellites is explained below.

## Identifying central subhaloes

Accretion of diffuse mass by the central subhalo happens before unbinding, meaning that we do not know *a priori* which subhalo in a halo is the most massive and hence which one is its central subhalo. Thus, HBT-HERONS needs to make an educated guess about which subhalo is the most likely to be the central one.

HBT-HERONS does this by identifying a mass threshold above which subhaloes are considered to be central candidates. The threshold is based on the most massive subhalo of the current halo in the previous output. HBT-HERONS uses the bound mass from the previous output because no unbinding has been performed yet, so masses at the current output are unknown. Thus, the threshold mass is:

$$
M_{\rm threshold} = f_{\rm major}M_{\rm 0}\,,
$$

where $f_{\rm major}$ is a user-defined parameter and $M_{\rm 0}$ the bound mass of the most massive subhalo of the current halo in the previous output. The effect of varying the value of $f_{\rm major}$ can be seen [here](./parameter_choices.md#majorprogenitormassratio).

Every subhalo whose bound mass in the previous output is lower than $M_{\rm threshold}$ is automatically classified to be a satellite. If more than one central subhalo candidate is found, HBT-HERONS chooses the candidate with the lowest orbital specific kinetic energy in the (current) frame of reference of the halo. Again, because unbinding has not been done yet, the centre of mass position and velocity of each central subhalo candidate uses the particles that were bound to it in the previous output.

## Emergent subhalo hierarchy

The process of identifying at each output time which subhaloes are centrals and which ones are satellites naturally leads to a subhalo hierarchy within a halo. Subhaloes that were centrals in the past had their own satellite systems, which remain connected to it even after they become a satellite subhalo. This means that connections reflect groups of subhaloes that fell in together because they were once their own independent halo.

Several [subhalo properties](../outputs/subhalo_properties.md#hierarchical-relationships) are saved by HBT-HERONS to keep track of the hierarchical relationships between subhaloes.