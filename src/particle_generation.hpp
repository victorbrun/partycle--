#ifndef PARTICLE_GENERATION_H
#define PARTICLE_GENERATION_H

#include "superellipsoid.hpp"
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
	// particle will be the same as reference_particle but with twice the volume.
	Superellipsoid reference_particle;
	Distribution volume_distribution;
	Distribution orientation_distribution;
};

std::vector<Superellipsoid*>* generate_random_particles(std::vector<ParticleDistribution> particle_distributions, 
														std::vector<double> volume_fractions, 
														int n_particles); 

Distribution parse_distribution(std::string distr_string);

/**
 * Computes the expected number of particles of each class 
 * needed to meet the target volume fraction. The returning
 * vector associates indeces to classes the same way as `target_volume_fractions`
 * and `expected_volumes`.
 *
 * @param `expected_volumes`: vector of expected value of volume for respective particle class.
 * @param `target_volume_fractions`: target volume fractions for each particle class.
 * @return the expected number of needed particle for each class to achive the target volume fractions
 */
std::vector<int> expected_particles_needed(std::vector<double> expected_volume_factors, std::vector<double> target_volume_fractions);

std::vector<double> expected_volumes(std::vector<ParticleDistribution> particle_distributions);

#endif
