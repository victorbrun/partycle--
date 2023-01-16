#include "utils.hpp"
#include "particle_generation.hpp"
#include "domain.hpp"

#include <eigen3/Eigen/src/SparseCore/SparsePermutation.h>
#include <iostream>
#include <random>

#define CONTACT_TOL 1e-2
#define COMPONENTS_FILE "../tests/example_particle_distributions.csv"

// Domain bounds 
#define X_LO 0
#define X_HI 10
#define Y_LO 0
#define Y_HI 10
#define Z_LO 0
#define Z_HI 10

int main() {
	
	std::vector<Component> components = parse_components(COMPONENTS_FILE);

	// Defines domain bounds 
	double x_range[2] = {X_LO, X_HI};
	double y_range[2] = {Y_LO, Y_HI};
	double z_range[2] = {Z_LO, Z_HI};
	double domain_vol = (x_range[1] - x_range[0]) * (y_range[1] - y_range[0]) * (z_range[1] - z_range[0]);

	// Generates particles acocrding to the specification given in component-file
	std::vector<Superellipsoid*>* particles = generate_random_particles(components, domain_vol);
	std::cout << "[INFO]: " << particles->size() << " particles generated" << std::endl;
	Domain domain = Domain(x_range, y_range, z_range, CONTACT_TOL, particles->size());

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> x(x_range[0], x_range[1]);
	std::uniform_real_distribution<double> y(y_range[0], y_range[1]);
	std::uniform_real_distribution<double> z(z_range[0], z_range[1]);
	

	// Randomly places sets centers of particles 
	for (size_t ix = 0; ix < particles->size(); ix++) {
		Eigen::Vector3d center(x(mt), y(mt), z(mt));	
		particles->at(ix)->set_center(center);
	}

	// TODO: Benchmark coordinate indexer.

	
}
