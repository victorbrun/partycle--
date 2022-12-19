#ifndef COORDINATE_INDEXER_H
#define COORDINATE_INDEXER_H

#include "superellipsoid.h"
#include <functional>

class CoordinateIndexer {
	private:
		Superellipsoid* particles;

		/**
		 * Indeces of particles in particles array sorted 
		 * by their x, y and z coordinates respectively.
		 */ 
		int* particles_idx_x_sorted;
		int* particles_idx_y_sorted;
		int* particles_idx_z_sorted;

		/**
		 * Performs binary seach on `a` under transformation by `func` to find the index at which `val` would
		 * be inserted to keep `a` ordered.
		 *
		 * @param `a`: array of ints ordered in ascending order defined by `func`, i.e. `func(a[i]) < func(a[i+1])`.
		 * @param `val`: double for which insert index is to be found. It is compared to elements in `a` after they have been transformed by `func`.
		 * @param `func`: funciton that defines the ordering of `a`.
		 * @return index where `val` can be inserted into `a` without fucking up the ordering of it.
		 */
		static int binary_search_insert_idx(int* a, double val, std::function<bool(double,double)> func);

		/**
		 * Performs binary search on `a` under transormation by `func` to find `val`.
		 *
		 * @param `a`: array of ints ordered in ascending order defined by `func`, i.e. `func(a[i]) < func(a[i+1])`.
		 * @param `val`: double to find in `a` after transormation under `func`.
		 * @param `func`: funciton that defines the ordering of `a`.
		 * @return index of `val` in `a`. If `val` is not present in `a` -1 is returned.
		 */
		static int binary_search(int* a, double val, std::function<bool(double,double)> func);
		
		/**
		 * Computes the intersection between two arrays of integers.
		 * The result is saved into `a1`. Note that if duplicate entries
		 * are present in either of the input arguments, unexpected behaviour
		 * may occur.
		 *
		 * @param 	`a1`, `a2` are arrays of integers.
		 */
		static void intersect(int* a1, int* a2);

	public:
};

#endif
