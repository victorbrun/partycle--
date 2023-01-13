#include "coordinate_indexer.hpp"
#include "superellipsoid.hpp"
#include <algorithm>
#include <cstdlib>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <unordered_set>

/***** Private *****/ 

int CoordinateIndexer::binary_search_insert_idx(std::vector<int>* a, double val, std::function<double(int)> func) {
	int hi = a->size()-1;
	int lo = 0;
	while (hi >= lo) {
		int mid = lo + (hi - lo)/2;
		if ( func(a->at(mid)) == val ) {
			return mid;
		} else if ( func(a->at(mid)) > val ) {
			hi = mid - 1;
		} else {
			lo = mid + 1;
		}
	}
	return hi + 1;
}

int CoordinateIndexer::binary_search(std::vector<int>* a, double val, std::function<double(int)> func) {return 0;}
		
std::vector<int> CoordinateIndexer::intersect(std::vector<int>& v1, std::vector<int>& v2) {
	if (v1.empty() || v2.empty()) return std::vector<int>();

	std::unordered_set<int> set {v1.begin(), v1.end()};
	std::vector<int> intersection;
	intersection.reserve(std::min(v1.size(), v2.size())); // Reserves largest possible space to increase speed a bit.
	for (int n: v2) {
		if (set.erase(n) > 0) { // if n exists in set, then 1 is returned and n is erased; otherwise, 0.
			intersection.push_back(n);
		}	
	}
	return intersection;
}

/***** Public *****/

CoordinateIndexer::CoordinateIndexer() {}

CoordinateIndexer::CoordinateIndexer(size_t n_particles) {
	this->particles = new std::vector<Superellipsoid*>;
	this->particles->reserve(n_particles);

	// Initiates the index vectors which will be n_particles large eventually
	this->particles_idx_x_sorted = new std::vector<int>;
	this->particles_idx_y_sorted = new std::vector<int>;
	this->particles_idx_z_sorted = new std::vector<int>;

	// Reserves space for vectors to make it faster to add entries to them.
	this->particles_idx_x_sorted->reserve(n_particles);
	this->particles_idx_y_sorted->reserve(n_particles);
	this->particles_idx_z_sorted->reserve(n_particles);
}

CoordinateIndexer::CoordinateIndexer(std::vector<Superellipsoid*>* particles) {
	this->particles = particles;

	// Initialise arrays to index the particles by their center coordniates
	size_t n_particles = particles->size();
	this->particles_idx_x_sorted = new std::vector<int>(n_particles);
	this->particles_idx_y_sorted = new std::vector<int>(n_particles);
	this->particles_idx_z_sorted = new std::vector<int>(n_particles);
	for (size_t ix = 0; ix < n_particles; ix++) {
		this->particles_idx_x_sorted->at(ix) = ix;
		this->particles_idx_y_sorted->at(ix) = ix;
		this->particles_idx_z_sorted->at(ix) = ix;
	}

	// Indexing the particles by their center coorindates
	std::sort(this->particles_idx_x_sorted->begin(), this->particles_idx_x_sorted->end(), [this](int i, int j) -> bool {
				float xi = this->particles->at(i)->get_center()[0];
				float xj = this->particles->at(j)->get_center()[0];
				return xi < xj;
			});	
	std::sort(this->particles_idx_y_sorted->begin(), this->particles_idx_y_sorted->end(), [this](int i, int j) -> bool {
				float yi = this->particles->at(i)->get_center()[1];
				float yj = this->particles->at(j)->get_center()[1];
				return yi < yj;
			});	
	std::sort(this->particles_idx_z_sorted->begin(), this->particles_idx_z_sorted->end(), [this](int i, int j) -> bool {
				float zi = this->particles->at(i)->get_center()[2];
				float zj = this->particles->at(j)->get_center()[2];
				return zi < zj;
			});	
}

CoordinateIndexer::~CoordinateIndexer() {
	std::cout << "[INFO]: destroying CoordinateIndexer instance" << std::endl;
	delete this->particles_idx_x_sorted;
	delete this->particles_idx_y_sorted;
	delete this->particles_idx_z_sorted;
}

void CoordinateIndexer::add_particle(Superellipsoid* p) {
	// The particles vector has no ordering so we just append
	this->particles->push_back(p);
	int new_idx = this->n_particles()-1;

	// Computes the place to insert the new particle 
	// to keep the indeces ordered
	Eigen::Vector3d c = p->get_center();
	int x_insert_idx = this->binary_search_insert_idx(this->particles_idx_x_sorted, c[0], [this](int ix) -> double {
				return this->get_particles()->at(ix)->get_center()[0];
			});
	int y_insert_idx = this->binary_search_insert_idx(this->particles_idx_y_sorted, c[1], [this](int ix) -> double {
				return this->get_particles()->at(ix)->get_center()[1];
			});
	int z_insert_idx = this->binary_search_insert_idx(this->particles_idx_z_sorted, c[2], [this](int ix) -> double {
				return this->get_particles()->at(ix)->get_center()[2];
			});

	// Inserts the new particle's index into each vector 
	particles_idx_x_sorted->insert(this->particles_idx_x_sorted->begin() + x_insert_idx, new_idx);
	particles_idx_y_sorted->insert(this->particles_idx_y_sorted->begin() + y_insert_idx, new_idx);
	particles_idx_z_sorted->insert(this->particles_idx_z_sorted->begin() + z_insert_idx, new_idx);
}

size_t CoordinateIndexer::n_particles() { 
	// length and index is of by one but next empty index is 
	// always one larger than last index to have element so this 
	// works great :)
	return this->particles->size();
}

std::vector<Superellipsoid*>* CoordinateIndexer::get_particles() {
	return this->particles;
}

std::vector<Superellipsoid*> CoordinateIndexer::particles_in_domain(double x_range[2], double y_range[2], double z_range[2]) {
	// Getting the indexes in the index vectors corresponding to the provided ranges.
	int x_lo = this->binary_search_insert_idx(this->particles_idx_x_sorted, x_range[0], [this](int ix) -> float {
				return this->particles->at(ix)->get_center()[0];				
			});
	int x_hi = this->binary_search_insert_idx(this->particles_idx_x_sorted, x_range[1], [this](int ix) -> float {
				return this->particles->at(ix)->get_center()[0];				
			});
	int y_lo = this->binary_search_insert_idx(this->particles_idx_x_sorted, y_range[0], [this](int ix) -> float {
				return this->particles->at(ix)->get_center()[1];				
			});
	int y_hi = this->binary_search_insert_idx(this->particles_idx_x_sorted, y_range[1], [this](int ix) -> float {
				return this->particles->at(ix)->get_center()[1];				
			});
	int z_lo = this->binary_search_insert_idx(this->particles_idx_x_sorted, z_range[0], [this](int ix) -> float {
				return this->particles->at(ix)->get_center()[2];				
			});
	int z_hi = this->binary_search_insert_idx(this->particles_idx_x_sorted, z_range[1], [this](int ix) -> float {
				return this->particles->at(ix)->get_center()[2];				
			});
	
	std::vector<int> x_particle_idx = std::vector<int>(this->particles_idx_x_sorted->begin() + x_lo, this->particles_idx_x_sorted->begin() + x_hi + 1);
	std::vector<int> y_particle_idx = std::vector<int>(this->particles_idx_y_sorted->begin() + y_lo, this->particles_idx_y_sorted->begin() + y_hi + 1);
	std::vector<int> z_particle_idx = std::vector<int>(this->particles_idx_z_sorted->begin() + z_lo, this->particles_idx_z_sorted->begin() + z_hi + 1);
	std::vector<int> xy_particle_idx = this->intersect(x_particle_idx, y_particle_idx);
	std::vector<int> xyz_particle_idx = this->intersect(xy_particle_idx, z_particle_idx);

	std::vector<Superellipsoid*> result(xyz_particle_idx.size());
	for (size_t ix = 0; ix < result.size(); ix++) {
		int idx = xyz_particle_idx.at(ix);
		result.at(ix) = this->particles->at(idx);	
	}

	return result;
}











