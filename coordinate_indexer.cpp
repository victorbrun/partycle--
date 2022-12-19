#include "coordinate_indexer.h"
#include "superellipsoid.h"

/***** Private *****/ 

int CoordinateIndexer::binary_search_insert_idx(int* a, double val, std::function<bool(double,double)> func) {return 0;}

int CoordinateIndexer::binary_search(int* a, double val, std::function<bool(double,double)> func) {return 0;}
		
void CoordinateIndexer::intersect(int* a1, int* a2) {}

/***** Public *****/

CoordinateIndexer::CoordinateIndexer() {}
CoordinateIndexer::CoordinateIndexer(int n_particles) {}
CoordinateIndexer::CoordinateIndexer(Superellipsoid* particles) {}

void CoordinateIndexer::insert(Superellipsoid p) {}

int CoordinateIndexer::n_particles() {return 0;}

Superellipsoid* CoordinateIndexer::get_particles() {
	// temp obv
	double scale[3] = {1,2,3};
	double shape[2] = {10,10};
	
	Superellipsoid p = Superellipsoid(0, scale, shape);
	return &p;
}

Superellipsoid* CoordinateIndexer::particles_in_domain(double x_range[2], double y_range[2], double z_range[2]) {
	return this->get_particles(); // temp obv	
}
