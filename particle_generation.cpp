#include "particle_generation.h"
#include "superellipsoid.h"
#include <algorithm>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <random>

std::vector<Superellipsoid*>* generate_random_particles(std::vector<ParticleDistribution> particle_distributions, 
														std::vector<double> volume_fractions, 
														int n_particles) 
{
	if ( particle_distributions.size() != volume_fractions.size() ) {
		throw std::invalid_argument("[ERROR]: size of particle_distributions does not equal size of volume_fractions.");
	}

	// Checks that each particle distribution has unique cls value since having two
	// distribution for the same class of particle is not well defined.
	for (int ix = 0; ix < particle_distributions.size(); ix++) {
		int cls1 = particle_distributions.at(ix).cls;
		for (int jx = ix+1; jx < particle_distributions.size(); ix++) {
			int cls2 = particle_distributions.at(jx).cls;
			if ( cls1 == cls2 ) {
				throw std::invalid_argument("[ERROR]: each ParticleDistribution must have a unique cls value.");
			}
		}
	}

	// Allocates vector on heap and generates random particles
	std::vector<Superellipsoid*>* particles = new std::vector<Superellipsoid*>(n_particles); 	
	for (int ix = 0; ix < n_particles; ix++) {
		
	}

	return particles;
}

std::vector<double> uniform_sampler(double lower_bound, double upper_bound) {
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(lower_bound, upper_bound);

	std::vector<double> v = {dist(mt)};
	return v;
}

std::vector<double> (*parse_distribution(std::string distr_string))(void) {
	int open_bracket_idx = distr_string.find("(");
	std::string distr_name = distr_string.substr(0, open_bracket_idx);
	std::string param_string = distr_string.substr(open_bracket_idx + 1, distr_string.size() - open_bracket_idx - 2);

	std::cout << "D = " << distr_name << ", params = " << param_string << std::endl;

	if (distr_name == "uniform") {
		int comma_idx = param_string.find(",");
		double lower_bound = std::stod(param_string.substr(0, comma_idx));
		double upper_bound = std::stod(param_string.substr(comma_idx + 1));
		return &uniform_sampler;
	} else if (distr_name == "normal") {

	} else if (distr_name == "log-normal") {

	} else {
		throw std::invalid_argument("[ERROR]: could not recognise distribution name. Available are: uniform, normal, log-normal.");
	}

}
