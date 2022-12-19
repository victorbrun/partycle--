#ifndef DOMAIN_H
#define DOMAIN_H

#include "superellipsoid.h"
#include "coordinate_indexer.h"

class Domain {
	private:
		

	public:
		// Constructor without particles or numbe rof particles to be added known
		Domain(float* x_range, float* y_range, float* z_range);
		// Constructor with particles available
		Domain(float* x_range, float* y_range, float* z_range, Superellipsoid* particles);
		// Constructor with number of particles to be added known
		Domain(float* x_range, float* y_range, float* z_range, int n_particles);

		// Getters
		
		// Setters

		void add_particle(Superellipsoid* p);
		Superellipsoid* particles_in_subdomain(float* x_range, float* y_range, float* z_range);
};

#endif
