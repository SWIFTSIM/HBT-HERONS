# Parameter choices

HBT-HERONS has several [free parameters](../running/parameter_file.md#subhaloes)
that can affect how subhaloes are tracked in time, subject to unbinding, and
which thresholds are used to determine if a subhalo is self-bound or not.

The default values that HBT-HERONS uses are intended to provide sensible results,
but the inquisitive user may wonder how changing each parameter would affect their
results. For this purpose, we provide in this page a series of tests where we
highlight how several relevant subhalo-based summary statistics vary when changing the value
of these key parameters.

The tests we provide were done using the intermediate resolution
($m_{\mathrm{dm}} = 6.72 \times 10^{9}\,\mathrm{M}_{\odot}$) DMO simulation of the $1\,\mathrm{Gpc}$ box
from the FLAMINGO suite of simulations ([Schaye et al 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.4978S/abstract)).

## Diagnostic plots

We show how a select number of summary statistics change when varying each free parameter that HBT-HERONS has. The summary statistics were chosen to quantify the instantaneous properties of the subhalo population, as well as their longer term time evolution. Our selection is by no means comprehensive, but it should provide the user the required intuition about how parameters influence the outcome of results.

The statistics are defined as follows:

* **Subhalo mass functions** ($M^{i}_{\mathrm{bound}}$): the number of self-bound subhaloes at $z = 0$, as a function of their current bound mass ($M^{\rm z=0}_{\mathrm{bound}}$) or their peak bound mass ($M^{\rm peak}_{\mathrm{bound}}$). The central *vs* satellite classification is done based on the `Depth` value of the subhalo at $z = 0$, where `Depth = 0` is a central subhalo and a satellite if `Depth != 0`.
* **Subhalo velocity functions** ($V^{i}_{\mathrm{max}}$): the number of self-bound subhaloes at $z = 0$, as a function of their current maximum circular velocity ($V^{\rm z=0}_{\mathrm{max}}$) or their peak circular velocity ($V^{\rm peak}_{\mathrm{max}}$). The central *vs* satellite classification is done based on the `Depth` value of the subhalo at $z = 0$, where `Depth = 0` is a central subhalo and a satellite if `Depth != 0`.
* **Length of the main progenitor branch** ($L_{\rm progenitor}$): how long the main progenitor branch of self-bound subhaloes at $z = 0$ is. The central *vs* satellite classification is done based on the `Depth` value of the subhalo at $z = 0$, where `Depth = 0` is a central subhalo and a satellite if `Depth != 0`. We further subdivide the subhalo population into three bins according to $M^{\rm z=0}_{\mathrm{bound}}$.
* **Normalised mass growth factor** ($\beta_{\rm M}$): the normalised arctangent of the exponent coefficient describing an exponential mass growth of a subhalo between two consecutive outputs, i.e.
$$
\beta_{\rm M} \equiv \dfrac{2}{\pi}\arctan\alpha_{\rm m} = \dfrac{2}{\pi}\arctan \dfrac{\ln [ M_{\rm bound}(z_{i+1}) / M_{\rm bound}(z_{i})]}{z_{i} - z_{i+1}}\,,
$$
where $z_{i} > z_{i+1}$. We use subhaloes across all redshifts, but only if both $M_{\rm bound}(z_{i+1})$ and $M_{\rm bound}(z_{i})$ are above a threshold mass $M_{\mathrm{th}} = 6.72\,\times10^{13}\,\mathrm{M}_{\odot}$ ($10^{4}$ particles). The central *vs* satellite classification is done based on the `Depth` value of the subhalo at the ***two consecutive outputs under consideration***. Thus, central subhaloes have `Depth = 0` in both outputs and satellite subhaloes have `Depth != 0` in both outputs instead.

Aside from these summary statistics, we also provide for certain parameters one or more supplementary plots, which are used to illustrate the effect that those parameters may have in a more direct manner.

!!! note

    We always keep the y-axis limits consistent across panels that show the same property, regardless of the parameter variation being shown. If less than 100 subhaloes are present in a given bin, we use a dashed line to indicate the poor statistics at that x-axis value.

<!-- These are left as future additions -->
<!-- #### `NumTracerHostFinding` -->
<!-- #### `NumTracersForDescendants` -->
<!-- #### `RefineMostBoundParticle` -->
<!-- #### `BoundFractionCenterRefinement` -->

#### `BoundMassPrecision`
=== "$M^{\rm z=0}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/mass_functions/z0/bound_mass_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/mass_functions/z0/bound_mass_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/mass_functions/z0/bound_mass_function_z0_satellites.png){ width="600" }
        </figure>

=== "$M^{\rm peak}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/mass_functions/peak/bound_mass_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/mass_functions/peak/bound_mass_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/mass_functions/peak/bound_mass_function_peak_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm z=0}_{\mathrm{max}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/velocity_functions/z0/vmax_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/velocity_functions/z0/vmax_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/velocity_functions/z0/vmax_function_z0_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm peak}_{\mathrm{max}}$"


    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/velocity_functions/peak/vmax_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/velocity_functions/peak/vmax_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/velocity_functions/peak/vmax_function_peak_satellites.png){ width="600" }
        </figure>

=== "$L_{\rm progenitor}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/main_branch_length/distribution_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/main_branch_length/distribution_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/main_branch_length/distribution_satellites.png){ width="600" }
        </figure>

=== "$\beta_{\rm M}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/growth_factors/growth_factor_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/growth_factors/growth_factor_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/BoundMassPrecision/growth_factors/growth_factor_satellites.png){ width="600" }
        </figure>

#### `SourceSubRelaxFactor`
=== "$M^{\rm z=0}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/mass_functions/z0/bound_mass_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/mass_functions/z0/bound_mass_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/mass_functions/z0/bound_mass_function_z0_satellites.png){ width="600" }
        </figure>

=== "$M^{\rm peak}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/mass_functions/peak/bound_mass_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/mass_functions/peak/bound_mass_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/mass_functions/peak/bound_mass_function_peak_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm z=0}_{\mathrm{max}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/velocity_functions/z0/vmax_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/velocity_functions/z0/vmax_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/velocity_functions/z0/vmax_function_z0_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm peak}_{\mathrm{max}}$"


    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/velocity_functions/peak/vmax_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/velocity_functions/peak/vmax_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/velocity_functions/peak/vmax_function_peak_satellites.png){ width="600" }
        </figure>

=== "$L_{\rm progenitor}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/main_branch_length/distribution_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/main_branch_length/distribution_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/main_branch_length/distribution_satellites.png){ width="600" }
        </figure>

=== "$\beta_{\rm M}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/growth_factors/growth_factor_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/growth_factors/growth_factor_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/SourceSubRelaxFactor/growth_factors/growth_factor_satellites.png){ width="600" }
        </figure>

#### `MaxSampleSizeOfPotentialEstimate`

=== "$M^{\rm z=0}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/mass_functions/z0/bound_mass_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/mass_functions/z0/bound_mass_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/mass_functions/z0/bound_mass_function_z0_satellites.png){ width="600" }
        </figure>

=== "$M^{\rm peak}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/mass_functions/peak/bound_mass_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/mass_functions/peak/bound_mass_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/mass_functions/peak/bound_mass_function_peak_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm z=0}_{\mathrm{max}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/velocity_functions/z0/vmax_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/velocity_functions/z0/vmax_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/velocity_functions/z0/vmax_function_z0_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm peak}_{\mathrm{max}}$"


    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/velocity_functions/peak/vmax_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/velocity_functions/peak/vmax_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/velocity_functions/peak/vmax_function_peak_satellites.png){ width="600" }
        </figure>

=== "$L_{\rm progenitor}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/main_branch_length/distribution_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/main_branch_length/distribution_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/main_branch_length/distribution_satellites.png){ width="600" }
        </figure>

=== "$\beta_{\rm M}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/growth_factors/growth_factor_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/growth_factors/growth_factor_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MaxSampleSizeOfPotentialEstimate/growth_factors/growth_factor_satellites.png){ width="600" }
        </figure>

#### `MinNumPartOfSub`

This parameter defaults to 20, and it is the minimum number of particles above which we consider a subhalo to be resolved. In our tests, we vary consider values of $N_{\mathrm{min}} \in \{10,20,32,50,100\}$.

The main difference is that the amplitude of the bound mass functions of resolved subhaloes at $z=0$ increases for subhaloes with $N_{\mathrm{bound}} > \mathrm{MinNumPartOfSub}$. This happens because HBT-HERONS uses an exclusive mass definition. (both at $z = 0$ and its peak value). In contrast,

This fact can easily be seen if we plot surface density of dark matter particles that are bound to a massive central subhalo of the box. For the largest value we use, $N_{\mathrm{min}} = 100$, several clumps. These correspond to dark matter. As we lower the threshold, these clumps disappear from the image because they are identified as self-bound, and hence their mass contribution the central is removed.

From the images, we also see that the clumps appear towards the outskirts, meaning that statistics such as the radial distribution functions will be similarly affected in a radial-dependent manner.
=== "$M^{\rm z=0}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/mass_functions/z0/bound_mass_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/mass_functions/z0/bound_mass_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/mass_functions/z0/bound_mass_function_z0_satellites.png){ width="600" }
        </figure>

=== "$M^{\rm peak}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/mass_functions/peak/bound_mass_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/mass_functions/peak/bound_mass_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/mass_functions/peak/bound_mass_function_peak_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm z=0}_{\mathrm{max}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/velocity_functions/z0/vmax_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/velocity_functions/z0/vmax_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/velocity_functions/z0/vmax_function_z0_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm peak}_{\mathrm{max}}$"


    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/velocity_functions/peak/vmax_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/velocity_functions/peak/vmax_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/velocity_functions/peak/vmax_function_peak_satellites.png){ width="600" }
        </figure>

=== "$L_{\rm progenitor}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/main_branch_length/distribution_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/main_branch_length/distribution_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/main_branch_length/distribution_satellites.png){ width="600" }
        </figure>

=== "$\beta_{\rm M}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/growth_factors/growth_factor_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/growth_factors/growth_factor_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinNumPartOfSub/growth_factors/growth_factor_satellites.png){ width="600" }
        </figure>


#### `MinCoreSize`


=== "$M^{\rm z=0}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/mass_functions/z0/bound_mass_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/mass_functions/z0/bound_mass_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/mass_functions/z0/bound_mass_function_z0_satellites.png){ width="600" }
        </figure>

=== "$M^{\rm peak}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/mass_functions/peak/bound_mass_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/mass_functions/peak/bound_mass_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/mass_functions/peak/bound_mass_function_peak_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm z=0}_{\mathrm{max}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/velocity_functions/z0/vmax_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/velocity_functions/z0/vmax_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/velocity_functions/z0/vmax_function_z0_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm peak}_{\mathrm{max}}$"


    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/velocity_functions/peak/vmax_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/velocity_functions/peak/vmax_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/velocity_functions/peak/vmax_function_peak_satellites.png){ width="600" }
        </figure>

=== "$L_{\rm progenitor}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/main_branch_length/distribution_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/main_branch_length/distribution_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/main_branch_length/distribution_satellites.png){ width="600" }
        </figure>

=== "$\beta_{\rm M}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/growth_factors/growth_factor_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/growth_factors/growth_factor_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MinCoreSize/growth_factors/growth_factor_satellites.png){ width="600" }
        </figure>

#### `MajorProgenitorMassRatio`

***Summary***: The instantaneous subhalo properties are not sensitive to the choice of `MajorProgenitorMassRatio`. However, this value can severely affect the time evolution of subhaloes. For values $\lesssim \mathcal{O}(0.1)$, the main progenitor branch length of well resolved subhaloes can be substantially shortened, due to low mass satellites being chosen as the central. Indeed, if one were to choose a value of 0, any subhalo (including orphans) could become the central. The distribution of mass growth factors are also extremely sensitive to changes in this value.


=== "$M^{\rm z=0}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/mass_functions/z0/bound_mass_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/mass_functions/z0/bound_mass_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/mass_functions/z0/bound_mass_function_z0_satellites.png){ width="600" }
        </figure>

=== "$M^{\rm peak}_{\mathrm{bound}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/mass_functions/peak/bound_mass_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/mass_functions/peak/bound_mass_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/mass_functions/peak/bound_mass_function_peak_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm z=0}_{\mathrm{max}}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/velocity_functions/z0/vmax_function_z0_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/velocity_functions/z0/vmax_function_z0_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/velocity_functions/z0/vmax_function_z0_satellites.png){ width="600" }
        </figure>

=== "$V^{\rm peak}_{\mathrm{max}}$"


    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/velocity_functions/peak/vmax_function_peak_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/velocity_functions/peak/vmax_function_peak_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/velocity_functions/peak/vmax_function_peak_satellites.png){ width="600" }
        </figure>

=== "$L_{\rm progenitor}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/main_branch_length/distribution_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/main_branch_length/distribution_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/main_branch_length/distribution_satellites.png){ width="600" }
        </figure>

=== "$\beta_{\rm M}$"

    === "All subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/growth_factors/growth_factor_all.png){ width="600" }
        </figure>

    === "Central subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/growth_factors/growth_factor_centrals.png){ width="600" }
        </figure>

    === "Satellite subhaloes"

        <figure markdown="span">
        ![image_title](../images/parameter_changes/MajorProgenitorMassRatio/growth_factors/growth_factor_satellites.png){ width="600" }
        </figure>