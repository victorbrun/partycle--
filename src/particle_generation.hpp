#ifndef PARTICLE_GENERATION_H
#define PARTICLE_GENERATION_H

#include "superellipsoid.hpp"
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <functional>
#include <random>
#include <vector>

struct Distribution {
	std::string name;
	std::vector<double> args;
};

struct Component {
	// Integer to differentiate between mixture components 
	int component_id;

	// Reference particle to which size distribution applies.
	// If a sample sample from volume_distribution yields 2, the resulting
	// particle will be the same as reference_particle but scaled to have volume 2.
	Superellipsoid reference_particle;
	Distribution volume_distribution;

	// In a mixture, this components should, on average, take up this much of the 
	// total volume of all the components (the volume of a component is the sum
	// of the volumes of every particle with the same class as a component).
	double target_volume_fraction;
};

/**
 * Generates enough random particles from `particle_distributions` such that 
 * their total volume is, on average, the same as domain volume. Furthermore, the fraction 
 * of the total volume of the particles, which is occupied by each particle class in 
 * `particle_distributions` is, on average, given by `volume_fractions`.
 *
 * Note that the generated particles have no initialised center, since this is assumed to 
 * to be set by the advanving front algoritm.
 *
 * @param `components`: vector specifying volume distribution, particle shape and volume fraction each component in mixture. The volume fraction is on average achived
 * @param `domain_volume`: the volume of the domain in which the generated particles is to be placed.
 */
std::vector<Superellipsoid*>* generate_random_particles(const std::vector<Component>& components, const double domain_volume); 

Distribution parse_distribution(const std::string& distr_string);

/**
 * Computes the expected number of particles of each component  
 * needed to meet the target volume fraction of respective component. The returning
 * vector is index aligned with `components`.
 *
 * @param `components`: vector of the mixtures different components.
 * @param `domain_volume`: volume of domain in which the particles will be places.
 * @return the expected number of needed particle for each component to achive the target volume fractions.
 */
std::vector<int> expected_particles_needed(const std::vector<Component>& components, const double domain_volume);

/**
 * Computes the probability with which to draw components from `components`
 * in order to on average have the volume fraction of each component in a large enough
 * set of particles be in accordance with their respective target volume fraction.
 *
 * @param `components`: vector of components.
 * @param `domain_volume`: volume of the domain in which the particles are to be placed.
 * @return vector of probabilities with which to select entries in `components` to achieve their respective target volume fractions.
 */
std::vector<double> particle_class_selection_prob(const std::vector<Component>& components,
												  const double domain_volume,
												  const size_t n_particles);

// Computes the analytical mean (expected value) of a random variable dsitributed according to `d`.
double mean(const Distribution& d);

/**
 * Initialises a sampler according to `d` and returns a function that samples from that distribution.
 *
 * @param `d`: distribution to create sampler from.
 * @param ´mt´: random number generator initialised with std::random_device.
 * @return parameter free function sampling from `d`.
 */
std::function<double(void)> get_sampler(const Distribution& d, std::mt19937& mt);

// Returns a quaternion where each rotation angle has been unifomrly sampled over its domain
Eigen::Quaternion<double> random_quaternion(std::mt19937& mt);

#endif
