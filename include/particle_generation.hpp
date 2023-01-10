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

struct ParticleDistribution {
	// Integer to differentiate between types of particles.
	int cls;

	// Reference particle to which size distribution applies.
	// If a sample sample from volume_distribution yields 2, the resulting
	// particle will be the same as reference_particle but scaled to have volume 2.
	Superellipsoid reference_particle;
	Distribution volume_distribution;
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
 * @param `particle_distributions`: vector specifying volume distribution and particle shape for each particle class.
 * @param `volume_fractions`: the volume fraction for each particle class in `particle_distributions`
 * 							  which is on average achived. This vector is index aligned with `particle_distributions`.
 * @param `domain_volume`: the volume of the domain in which the generated particles is to be placed.
 */
std::vector<Superellipsoid*>* generate_random_particles(const std::vector<ParticleDistribution>& particle_distributions, 
														const std::vector<double>& volume_fractions, 
														const double domain_volume); 

Distribution parse_distribution(std::string distr_string);

/**
 * Computes the expected number of particles of each class 
 * needed to meet the target volume fraction. The returning
 * vector associates indeces to classes the same way as `target_volume_fractions`
 * and `expected_volumes`.
 *
 * @param `expected_volumes`: vector of expected value of volume for respective particle class.
 * @param `target_volume_fractions`: target volume fractions for each particle class.
 * @param `domain_volume`: volume of domain in which the particles will be places.
 * @return the expected number of needed particle for each class to achive the target volume fractions
 */
std::vector<int> expected_particles_needed(const std::vector<ParticleDistribution>& particle_distributions, 
										   const std::vector<double>& target_volume_fractions,
										   const double domain_volume);

/**
 * Computes the probability with which to draw particle classes from `particle_distributions`
 * in order to on average have the volume fraction of each particle class in a large enough
 * set of particles be in accordance with `target_volume_fractions`.
 *
 * @param `particle_distributions`: vector of distributions for the different particle classes.
 * @param `target_volume_fractions`: the target volume fraction for each particles class sought to achieve.
 * 									 Index aligned with `particle_distributions`.
 * @param `domain_volume`: volume of the domain in which the particles are to be placed.
 * @return vector of probabilities with which to select entries in `particle_distributions`to achieve `target_volume_fractions`.
 */
std::vector<double> particle_class_selection_prob(const std::vector<ParticleDistribution>& particle_distributions,
												  const std::vector<double>& target_volume_fractions,
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

// Returns a quaternion where each rotation angle has been unifomrly sampled  over its domain
Eigen::Quaternion<double> random_quaternion(std::mt19937& mt);

#endif
