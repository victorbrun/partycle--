#ifndef COORDINATE_INDEXER_H
#define COORDINATE_INDEXER_H

#include "superellipsoid.hpp"
#include <functional>
#include <vector>

class CoordinateIndexer {
	private:
		std::vector<Superellipsoid*>* particles;
		/**
		 * Indeces of particles in particles array sorted 
		 * by their x, y and z coordinates respectively.
		 */ 
		std::vector<int>* particles_idx_x_sorted;
		std::vector<int>* particles_idx_y_sorted;
		std::vector<int>* particles_idx_z_sorted;

		/**
		 * Performs binary seach on `a` under transformation by `func` to find the index at which `val` would
		 * be inserted to keep `a` ordered.
		 *
		 * @param `a`: pointer to vector of ints ordered in ascending order defined by `func`, i.e. `func(a[i]) < func(a[i+1])`.
		 * @param `val`: double for which insert index is to be found. It is compared to elements in `a` after they have been transformed by `func`.
		 * @param `func`: funciton that defines the ordering of `a`.
		 * @return index where `val` can be inserted into `a` without fucking up the ordering of it.
		 */
		static int binary_search_insert_idx(std::vector<int>* a, double val, std::function<double(int)> func);

		/**
		 * Performs binary search on `a` under transormation by `func` to find `val`.
		 *
		 * @param `a`: pointer to vector of ints ordered in ascending order defined by `func`, i.e. `func(a[i]) < func(a[i+1])`.
		 * @param `val`: double to find in `a` after transormation under `func`.
		 * @param `func`: funciton that defines the ordering of `a`.
		 * @return index of `val` in `a`. If `val` is not present in `a` -1 is returned.
		 */
		static int binary_search(std::vector<int>* a, double val, std::function<double(int)> func);
		
		/**
		 * Computes the intersection between two arrays of integers.
		 * The result is saved into `a1`. Note that if duplicate entries
		 * are present in either of the input arguments, unexpected behaviour
		 * may occur.
		 *
		 * @param `a1`, `a2` references to vectors of integers.
		 */
		static std::vector<int> intersect(std::vector<int>& a1, std::vector<int>& a2);

	public:
		// Constructor for when neither number of particles or particles are known.
		CoordinateIndexer();
		
		// Constructor for when particles are known.
		CoordinateIndexer(std::vector<Superellipsoid*>* particles);
		
		// Constructor for when number of particles are known.
		CoordinateIndexer(int n_particles);

		// Frees all the memory allocated on the heap
		void destroy();

		// Adds `p` to CoordinateIndexer and indexes its coordinates.
		void add_particle(Superellipsoid* p);

		// Returns array of all particles.
		std::vector<Superellipsoid*>* get_particles();
	
		// Returns number of particles in CoordinateIndexer.
		int n_particles();

		/**
		 * Returns array of particles whose center lie in the domain defined by `x_range`x`y_range`x`z_range`.
		 *
		 * @param `x_range`: array of 2 doubles defining x-axis range. Smallest value is at index 0 and larges is at index 1.
		 * @param `y_range`: array of 2 doubles defining y-axis range. Smallest value is at index 0 and larges is at index 1.
		 * @param `z_range`: array of 2 doubles defining z-axis range. Smallest value is at index 0 and larges is at index 1.
		 * @return vector of pointers to superellipsoids all of which lie in the domain defined by `x_range`x`y_range`x`z_range`.
		 */
		std::vector<Superellipsoid*> particles_in_domain(double x_range[2], double y_range[2], double z_range[2]);
};

#endif
