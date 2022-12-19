#include "superellipsoid.h"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <iostream>
#include <ostream>

int main() {
	float scale_params[3] = {1,1,1};
	float shape_params[2] = {2,2};
	Superellipsoid p = Superellipsoid(0, scale_params, shape_params);
	Eigen::Vector3d center(2,2,2);
	Eigen::Quaternion<double> rot(1,1,1,1);
	
	p.set_center(center);
	p.set_orientation(rot);


	Eigen::Vector3d x(10,10,10);
	std::cout<<"INSIDE-OUTSIDE: "<<p.inside_outside_grad(x)<<std::endl;

	return 0;
}
