# Unbinding and mass stripping


### Particle subsampling

To speed up the process of unbinding, HBT-HERONS restricts the number of particles that are gravity sources to a user-defined maximum quantity ($N_{\rm subsample}$; [`MaxSampleSizeOfPotentialEstimate`](../running/parameter_file.md#unbinding)). If the number of particles in the source subhalo is greater than the chosen threshold, a subsample of $N_{\rm subsample}$ particles is chosen at random. Only the chosen particles contribute to the gravitational potential energy, with all other particles acting as test particles of the gravitational potential field. The mass of the randomly chosen subset of individual particles is scaled to ensure mass conservation.

Note that there are additional options that can refine the way in which the random subsampling occurs. For example, the subsampling of specific particle types (e.g. black holes) can be explicitly disabled. Moreover, the scaling of particle masses can be done on a per-particle-type basis, which is recommended if particle types are spatial distributed in very different ways with respect to each another.

We show below how the subsampling of particles works in practice, and how it may bias results compared to not subsampling, using an idealised subhalo.

<h4>Isolated subhalo example</h4>

We generate a spherically symmetric subhalo in isolation and equilibrium using $N = 10^{5}$ particles. Since the number of particles is above the default threshold of $N_{\rm subsample} = 10^{3}$, we randomly select $10^{3}$ particles. The projection of the true (blue dots) and subsampled (red dots) spatial distributions are shown below.

<figure markdown="span">
![image_title](../images/algorithm/idealised_subhalo_sampling.png){ width="600" }
</figure>

We can visually see that the distribution of randomly chosen subset of particles *approximately* follows the underlying ground truth. The extent is however different, and regions with few particles become prone to shot noise as we reduce the number of points representing the mass distribution. These problems exacerbate as the chosen value of $N_{\rm subsample}$ decreases, since you become increasingly reliant on whether your random choice of particles is representative or not.

However, when we perform the subsampling, we are only interested in speeding up the process of identifying *which particles are self-bound*. Hence, the more appropiate metric of how well subsampling works is to investigate how the number of bound particles <br> ($N_{\rm bound}$) differs from the ground truth. We show the cumulative distribution of $N_{\rm bound}$ for 1000 different random realisations for $N_{\rm subsample} = 10, 10^{2}$ and $10^{3}$. Since this is an isolated subhalo in equilibrium, all particles should be bound ($N_{\rm bound} = 10^{5}$).

<figure markdown="span">
![image_title](../images/algorithm/subsampling_distribution_Nbound.png){ width="600" }
</figure>

For this test, which represents a best case scenario because of the spherical symmetry, we see that using as few as $10^{2}$ particles yields results comparable to using all particles. The only difference is that the amount of calculations has drastically reduced, and hence the computational cost is lower. However, using very few particles ($N_{\rm subsample} = 10$) results in $\approx 8\%$ of cases in an underestimate of the bound mass of the subhalo. The cases where the bound mass is underestimated reflect a random subsample that poorly represents the underlying mass distribution, which becomes more likely to happen as the ratio $N_{\rm bound} / N_{\rm subsample}$ becomes very small.
