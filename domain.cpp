#include "domain.h"
#include "superellipsoid.h"
#include <algorithm>

/***** Private *****/

void Domain::remove_particle_af(Superellipsoid* p) { 
	// Searches for the index at which p is located in advancing_front
	int remove_idx = -1;
	for (int ix = 0; ix < this->advancing_front.size(); ix++) {
		if (this->advancing_front.at(ix) == p) remove_idx = ix;
	}

	// Only remove element if p was found in af
	if (remove_idx > -1) this->advancing_front.erase(this->advancing_front.begin() + remove_idx); 
}

void Domain::add_particle_af(Superellipsoid* p) {
	this->advancing_front.push_back(p);
}

std::vector<Superellipsoid*> Domain::particles_in_subdomain(double x_range[2], double y_range[2], double z_range[2]) {
	return this->particles.particles_in_domain(x_range, y_range, z_range);
}

int binary_approach(Superellipsoid* fixed_particles[3], Superellipsoid* mobile_particle) { return 0; }

/***** Public *****/

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[1];
	this->z_bounds[1] = z_range[1];

	this->particles = CoordinateIndexer();
}

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, int n_particles) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[1];
	this->z_bounds[1] = z_range[1];

	this->particles = CoordinateIndexer(n_particles);
}

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, std::vector<Superellipsoid*>* particles) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[1];
	this->z_bounds[1] = z_range[1];

	this->particles = CoordinateIndexer(particles);
}

int Domain::n_particles() { return this->particles.n_particles(); }
