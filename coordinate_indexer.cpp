#include "coordinate_indexer.h"
#include "superellipsoid.h"
#include <algorithm>
#include <cstdlib>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <iostream>

/***** Private *****/ 

int CoordinateIndexer::binary_search_insert_idx(int* a, double val, std::function<bool(double,double)> func) {return 0;}

int CoordinateIndexer::binary_search(int* a, double val, std::function<bool(double,double)> func) {return 0;}
		
void CoordinateIndexer::intersect(int* a1, int* a2) {}

/***** Public *****/

CoordinateIndexer::CoordinateIndexer() {}

CoordinateIndexer::CoordinateIndexer(int n_particles) {}

CoordinateIndexer::CoordinateIndexer(Superellipsoid** particles, int n_particles) {
	this->particles = particles;
	this->particles_len = n_particles;
	this->next_empty_index = this->n_particles();

	// Initialise arrays to index the particles by their center coordniates
	this->particles_idx_x_sorted = new std::vector<int>(n_particles);
	this->particles_idx_y_sorted = new std::vector<int>(n_particles);
	this->particles_idx_z_sorted = new std::vector<int>(n_particles);
	for (int ix = 0; ix < n_particles; ix++) {
		this->particles_idx_x_sorted->at(ix) = ix;
		this->particles_idx_y_sorted->at(ix) = ix;
		this->particles_idx_z_sorted->at(ix) = ix;
	}

	// Indexing the particles by their center coorindates
	std::sort(this->particles_idx_x_sorted->begin(), this->particles_idx_x_sorted->end(), [this](int i, int j){
				float xi = this->particles[i]->get_center()[0];
				float xj = this->particles[j]->get_center()[0];
				return xi < xj;
			});	
	std::sort(this->particles_idx_y_sorted->begin(), this->particles_idx_y_sorted->end(), [this](int i, int j){
				float yi = this->particles[i]->get_center()[1];
				float yj = this->particles[j]->get_center()[1];
				return yi < yj;
			});	
	std::sort(this->particles_idx_z_sorted->begin(), this->particles_idx_z_sorted->end(), [this](int i, int j){
				float zi = this->particles[i]->get_center()[2];
				float zj = this->particles[j]->get_center()[2];
				return zi < zj;
			});	

	// Printing for debugg
	for (int ix = 0; ix < n_particles; ix++) {
		std::cout<<particles_idx_x_sorted->at(ix)<<" ";
	}
	std::cout<<"\n";
}

void CoordinateIndexer::insert(Superellipsoid* p) {
	if (!this->can_insert) throw std::domain_error("ERROR: cannot insert new particle when object have been initialised using current constructor");

	// Adds particle to particles array and does some houskeeping
	// to keep everything up to date
	int new_index = this->next_empty_index;
	this->particles[new_index] = p;
	this->particles_len++;
	this->next_empty_index++;

	// Extarcts the coordinates of new particles 
	Eigen::Vector3d c = p->get_center();
	double x = c[0];
	double y = c[1];
	double z = c[2];

	
}

void CoordinateIndexer::destroy() {
	delete [] this->particles;
	delete [] this->particles_idx_x_sorted;
	delete [] this->particles_idx_y_sorted;
	delete [] this->particles_idx_z_sorted;
}


int CoordinateIndexer::n_particles() { 
	// length and index is of by one but next empty index is 
	// always one larger than last index to have element so this 
	// works great :)
	return this->next_empty_index;
}

Superellipsoid** CoordinateIndexer::get_particles() {
	return this->particles;
}

Superellipsoid** CoordinateIndexer::particles_in_domain(double x_range[2], double y_range[2], double z_range[2]) {
	return this->get_particles(); // temp obv	
}
