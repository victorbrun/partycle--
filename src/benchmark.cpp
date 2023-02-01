#include "superellipsoid.hpp"
#include "utils.hpp"
#include "particle_generation.hpp"
#include "domain.hpp"
#include "particle_generation.hpp"

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

// Benchmark settings
constexpr size_t n_samples = 10;
constexpr size_t n_iter = 1000;

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

BASELINE(increment_front, baseline, n_samples, n_iter) {
	ref_p1_scale[3] = {1, 1, 1};
	ref_p1_shape[2] = {2, 2};
	ref_p1 = Superellipsoid(1, ref_p1_scale, ref_p1_shape);

	std::vector<Component> components = {
		{1, },
		
	}

	// Defines domain bounds 
	double x_range[2] = {X_LO, X_HI};
	double y_range[2] = {Y_LO, Y_HI};
	double z_range[2] = {Z_LO, Z_HI};
	double domain_vol = (x_range[1] - x_range[0]) * (y_range[1] - y_range[0]) * (z_range[1] - z_range[0]);

	// Generates particles acocrding to the specification given in component-file
	std::vector<Superellipsoid*>* particles = generate_random_particles(components, domain_vol);

	// Shuffles the particles so that by iterating over the vector 
	// will be the same as randomly selecting a component and generating 
	// a particle according to it. 
	// TODO: check that this argument of shuffling pre-generated particles is equiv. to that of randomly sampling
	// component and then generating particle.
	std::random_device rd;
	std::mt19937 mt(rd());
	std::shuffle(particles->begin(), particles->end(), mt);

	// Initialises domain. This is done here and not where its bounds are defined since we want to 
	// specify the number of particles we are going to add to it. This will make the domain pre-allocate
	// memory.
	Domain domain = Domain(x_range, y_range, z_range, contact_tol, particles->size());
	
	// Uses first four particles to initiate advancing front
	Superellipsoid* init_p[4] = {particles->at(0), particles->at(1), particles->at(2), particles->at(3)};
	domain.initialise_outward_advancing_front(init_p);
	
	// Uses the rest of the particles to incerement the advanding front
	for (size_t ix = 4; ix < particles->size(); ix++) {
		domain.increment_advancing_front(particles->at(ix));
	}
}






