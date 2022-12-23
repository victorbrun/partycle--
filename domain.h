#ifndef DOMAIN_H
#define DOMAIN_H

#include "superellipsoid.h"
#include "coordinate_indexer.h"
#include <eigen3/Eigen/Sparse>
#include <vector>

class Domain {
	private:
		// Class that contains all particles
		CoordinateIndexer particles;
		
		// Bounds of the domain
		double x_bounds[2];
		double y_bounds[2];
		double z_bounds[2];

		// Advancing front
		std::vector<Superellipsoid*> advancing_front;

		// Neighbouring particles stuff
		double contact_tol; // if |distance(p1,p2)| <= contact_tol p1 and p2 are considered to be in contact 
		Eigen::SparseMatrix<int> contact_matrix; 

		// Removes removes the given particle pointer from the advancing front.
		void remove_particle_af(Superellipsoid* p);
	
		/**
		 * Performes binary approach moving `mobile_particle` toward the average of the centers of `fixed_particles`.
		 * The function returns an int specifying the exit state.
		 *
		 * @param `fixed_particles`: array of three pointers to Superellipsoids.
		 * @param `mobile_particle`: pointer to particle which is moved.
		 * @return 
		 */
		int binary_approach(Superellipsoid* fixed_particles[3], Superellipsoid* mobile_particle);

		std::vector<Superellipsoid*> particles_in_subdomain(float* x_range, float* y_range, float* z_range);
	
	public:
		// Constructor without particles or numbe rof particles to be added known
		Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol);
		// Constructor with particles available
		Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, std::vector<Superellipsoid*>* particles);
		// Constructor with number of particles to be added known
		Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, int n_particles);

		int n_particles();
		
		// Adds `p` to the domain
		void add_particle(Superellipsoid* p);
		
		void initialise_advancing_front(std::string direction, std::vector<Superellipsoid*> particles);

		// Introduces the given particle into the domain and increments the advancing front accordingly.
		void increment_advancing_front(Superellipsoid* p);

		void draw(int samples);

};

#endif
