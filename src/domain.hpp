#ifndef DOMAIN_H
#define DOMAIN_H

#include "superellipsoid.hpp"
#include "coordinate_indexer.hpp"
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

		// Up to date record of sphere radius which is guranteed to enclose all particles in domain 
		double larges_circumscribing_sphere_radius;

		// Removes removes the given particle pointer from the advancing front.
		void remove_particle_af(Superellipsoid* p);

		/**
		 * Adds particle to advancing front. It does not check
		 * if `p` is already in advancing_front, it can thus
		 * produce an advancing_front which has several entries
		 * pointing to the same superellipsoid.
		 *
		 * @param `p`: pointer to the superellipsoid which is to be added to advancing_front.
		 */
		void add_particle_af(Superellipsoid* p);

		/**
		 * Initialises the an outward going advancing front by placing the 
		 * four particles in `particles` in the corners of a tetrahedron with center in 
		 * the domains midpoint. The superellipsoids in `particles` may not be present
		 * in domain or advancing_front since they will be added there-too by this function.
		 *
		 * @param `particles`: four pointers to superellipsoids. Center of the superellipsoids will be changed.
		 */
		void initialise_outward_advancing_front(Superellipsoid* particles[4]);

		/**
		 * Performes binary approach moving `mobile_particle` as close to the average of the centers of `fixed_particles`
		 * without overlapping any particle. The function returns an int specifying the exit state.
		 *
		 * @param `fixed_particles`: array of three pointers to Superellipsoids.
		 * @param `mobile_particle`: pointer to particle which is moved.
		 * @return 
		 * 		0: success: distance between mobile particle and any fixed particle < tol.
		 * 		1: --
		 * 		2: --
		 * 		3: --
		 */
		int binary_approach(Superellipsoid* fixed_particles[3], Superellipsoid* mobile_particle);

		/**
		 * Returns array of particles whose center lie in the domain defined by `x_range`x`y_range`x`z_range`.
		 *
		 * @param `x_range`: array of 2 doubles defining x-axis range. Smallest value is at index 0 and larges is at index 1.
		 * @param `y_range`: array of 2 doubles defining y-axis range. Smallest value is at index 0 and larges is at index 1.
		 * @param `z_range`: array of 2 doubles defining z-axis range. Smallest value is at index 0 and larges is at index 1.
		 * @return vector of pointers to superellipsoids all of which lie in the domain defined by `x_range`x`y_range`x`z_range`.
		 */
		std::vector<Superellipsoid*> particles_in_subdomain(double* x_range, double* y_range, double* z_range);
	
	public:
		// Constructor without particles or numbe rof particles to be added known
		Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol);
		
		// Constructor with number of particles to be added known
		Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, int n_particles);

		// Constructor with particles available
		Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, std::vector<Superellipsoid*>* particles);
		
		/** 
		 * Reserves space for n_particles. Ought to be used when you need domain
		 * volume in order to compute the number of particles if will be filled with.
		 */
		void reserve_particle_space(int n_particles);

		int n_particles(void);
		
		// Adds `p` to the domain
		void add_particle(Superellipsoid* p);
		
		void initialise_advancing_front(std::string direction, std::vector<Superellipsoid*> particles);

		// Introduces the given particle into the domain and increments the advancing front accordingly.
		void increment_advancing_front(Superellipsoid* p);

		// Returns the volume of the domain
		double volume(void);

		void draw(int samples);

};

#endif
