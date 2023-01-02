#include "superellipsoid.h"
#include "coordinate_indexer.h"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
	int n_particles = 100;
	std::vector<Superellipsoid*>* super_arr = new std::vector<Superellipsoid*>(n_particles);
	for (int ix = 0; ix < n_particles; ix++) {
		double scale_params[3] = {1,1,1};
		double shape_params[2] = {10,20};
		super_arr->at(ix) = new Superellipsoid(0, scale_params, shape_params);
		Eigen::Vector3d center(100-ix,100-ix,100-ix);
		Eigen::Quaternion<double> rot(1,1,1,1);

		super_arr->at(ix)->set_center(center);
		super_arr->at(ix)->set_orientation(rot);
	}

	CoordinateIndexer ci = CoordinateIndexer(super_arr);	
	
	double x_range[2] = {4, 20};
	double y_range[2] = {4, 20};
	double z_range[2] = {4, 20};
	std::vector<Superellipsoid*> result = ci.particles_in_domain(x_range, y_range, z_range);
	
	for (int ix = 0; ix < result.size(); ix++) {
		auto c = result[ix]->get_center();
		std::cout << "ix =" << ix << ": (x,y,z) = (" << c[0] << "," << c[1] << "," << c[2] << ")" << std::endl;
	}

	for (int ix = 0; ix < n_particles; ix++) {
		double scale_params[3] = {1,1,1};
		double shape_params[2] = {10,20};
		Eigen::Vector3d center(100-ix,100-ix,100-ix);
		Eigen::Quaternion<double> rot(1,1,1,1);

		Superellipsoid* p = new Superellipsoid(0, scale_params, shape_params);
		p->set_center(center);
		p->set_orientation(rot);
		
		ci.add_particle(p);
	}

	std::cout << "\n";
	result = ci.particles_in_domain(x_range, y_range, z_range);	
	for (int ix = 0; ix < result.size(); ix++) {
		auto c = result[ix]->get_center();
		std::cout << "ix =" << ix << ": (x,y,z) = (" << c[0] << "," << c[1] << "," << c[2] << ")" << std::endl;
	}

	return 0;
}
