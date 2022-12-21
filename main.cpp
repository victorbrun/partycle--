#include "superellipsoid.h"
#include "coordinate_indexer.h"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <iostream>

int main() {
	int n_particles = 10;
	Superellipsoid* super_arr[10];
	for (int ix = 0; ix < n_particles; ix++) {
		double scale_params[3] = {1,1,1};
		double shape_params[2] = {10,20};
		super_arr[ix] = new Superellipsoid(0, scale_params, shape_params);
		Eigen::Vector3d center(10-ix,0,0);
		Eigen::Quaternion<double> rot(1,1,1,1);

		super_arr[ix]->set_center(center);
		super_arr[ix]->set_orientation(rot);
	}

	CoordinateIndexer ci = CoordinateIndexer(super_arr, n_particles);	


	return 0;
}
