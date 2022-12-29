#ifndef PARTICLE_GENERATION_H
#define PARTICLE_GENERATION_H

#include "superellipsoid.h"
#include <vector>

struct ParticleDistribution {
	int cls;
	std::string size_distribution;
	std::string orientation_distribution;
};

std::vector<Superellipsoid*>* generate_random_particles(std::vector<ParticleDistribution> particle_distributions, 
														std::vector<double> volume_fractions, 
														int n_particles); 

std::vector<double> (*parse_distribution(std::string distr_string))(void);


#endif
