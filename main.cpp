#include "superellipsoid.h"
#include "coordinate_indexer.h"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <iostream>
#include <vector>

int main() {
	int n_particles = 10;
	Superellipsoid* super_arr[10];
	for (int ix = 0; ix < n_particles; ix++) {
		double scale_params[3] = {1,1,1};
		double shape_params[2] = {10,20};
		super_arr[ix] = new Superellipsoid(0, scale_params, shape_params);
		Eigen::Vector3d center(ix,ix,ix);
		Eigen::Quaternion<double> rot(1,1,1,1);

		super_arr[ix]->set_center(center);
		super_arr[ix]->set_orientation(rot);
	}

	CoordinateIndexer ci = CoordinateIndexer(super_arr, n_particles);	
	
	double x_range[2] = {4, 6};
	double y_range[2] = {4, 6};
	double z_range[2] = {4, 6};
	std::vector<Superellipsoid*> result = ci.particles_in_domain(x_range, y_range, z_range);
	
	
	for (int ix = 0; ix < result.size(); ix++) {
		auto c = result[ix]->get_center();
		std::cout << "ix =" << ix << ": (x,y,z) = (" << c[0] << "," << c[1] << "," << c[2] << ")" << std::endl;
	}

	return 0;
}
