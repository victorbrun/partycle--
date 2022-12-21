#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include "superellipsoid.h"

// Sign function. Returns the sign of input and 0 if input is 0.
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

Superellipsoid::Superellipsoid(int cls, double scale_params[3], double shape_params[2]) {
	double n1 = shape_params[0];
	double n2 = shape_params[1];
	if (n1 < 2 || n2 < 2) {
		throw std::invalid_argument("ERROR: n1 and n2 must be larger or equal to 2 to ensure that superellipsoid surface is two times differentiable");
	}
	this->shape[0] = n1;
	this->shape[1] = n2;

	double a = scale_params[0];
	double b = scale_params[1];
	double c = scale_params[2];
	this->scale[0] = a;
	this->scale[1] = b;
	this->scale[2] = c;
		
	this->cls = cls;
}

int Superellipsoid::get_class() { return this->cls; }
Eigen::Vector3d Superellipsoid::get_center() { return this->c; }
Eigen::Quaternion<double> Superellipsoid::get_orientation() { return this->q; }

double Superellipsoid::get_scale(std::string param_name) {
	if (param_name == "a") {
		return scale[0];
	} else if (param_name == "b") {
		return scale[1];
	} else if (param_name == "c") {
		return scale[2];
	} else {
		throw std::invalid_argument("ERROR: get_scale accepts param_name: \"a\", \"b\", \"c\"");
	}
}

double Superellipsoid::get_shape(std::string param_name) {
	if (param_name == "n1") {
		return shape[0];
	} else if (param_name == "n2") {
		return shape[1];
	} else {
		throw std::invalid_argument("ERROR: get_shape accepts param_name: \"n1\", \"n1\"");
	}
}

void Superellipsoid::set_center(Eigen::Vector3d c) { this->c = c; }
void Superellipsoid::set_orientation(Eigen::Quaternion<double> q) { this->q = q; }

Eigen::Vector3d Superellipsoid::to_local_coords(Eigen::Vector3d x) {
	Eigen::Matrix3d R = this->get_orientation().toRotationMatrix();
	return  R * ( x - this->get_center() );
}

double Superellipsoid::implicit_surface(Eigen::Vector3d x) {
	Eigen::Vector3d x_local = this->to_local_coords(x);
	double xl = x_local[0];
	double yl = x_local[1];
	double zl = x_local[2];

	double a 	= this->get_scale("a"); 
	double b 	= this->get_scale("b"); 
	double c 	= this->get_scale("c");
	double n1 	= this->get_shape("n1");
	double n2 	= this->get_shape("n2");

	return std::pow( std::pow(std::abs(xl/a), n2) + std::pow(std::abs(yl/b), n2), n1/n2 ) + std::pow(std::abs(zl/c), n1);

	return 0;
}

double Superellipsoid::inside_outside(Eigen::Vector3d x) {
	return this->implicit_surface(x) - 1;
}

Eigen::Vector3d Superellipsoid::inside_outside_grad(Eigen::Vector3d x) {
	Eigen::Vector3d x_local = this->to_local_coords(x);
	double xl = x_local[0];
	double yl = x_local[1];
	double zl = x_local[2];

	double a 	= this->get_scale("a"); 
	double b 	= this->get_scale("b"); 
	double c 	= this->get_scale("c");
	double n1 	= this->get_shape("n1");
	double n2 	= this->get_shape("n2");

	double v = std::pow(std::abs(xl/a), n2) + std::pow(std::abs(yl/b), n2);
	double dx = (n1/a) * std::pow(std::abs(xl/a), n2 - 1) * std::pow(v, n1/n2 - 1) * sgn(xl);
	double dy = (n1/b) * std::pow(std::abs(yl/b), n2 - 1) * std::pow(v, n1/n2 - 1) * sgn(yl);
	double dz = (n1/c) * std::pow(std::abs(zl/c), n1 - 1) * sgn(zl);

	Eigen::Vector3d grad(dx, dy, dz);
	return grad;
}

Eigen::Matrix3d Superellipsoid::inside_outside_hess(Eigen::Vector3d x) {
	Eigen::Vector3d x_local = this->to_local_coords(x);
	double xl = x_local[0];
	double yl = x_local[1];
	double zl = x_local[2];

	double a 	= this->get_scale("a"); 
	double b 	= this->get_scale("b"); 
	double c 	= this->get_scale("c");
	double n1 	= this->get_shape("n1");
	double n2 	= this->get_shape("n2");

	double v = std::pow(std::abs(xl/a), n2) + std::pow(std::abs(yl/b), n2);
	
	double dxdx = (n1*(n2-1)/(a*a)) * std::pow(std::abs(xl/a), n2 - 2) * std::pow(v, n1/n2 - 1) + (n2*(n1-n2)/(a*a)) * std::pow(std::abs(xl/a), 2*n2 - 2) * std::pow(v, n1/n2 - 2);
	double dxdy = (n1*(n1-n2)/(a*b)) * std::pow(std::abs(xl/a), n2 - 1) * std::pow(std::abs(yl/b), n2 - 1) * std::pow(v, n1/n2 - 2) * sgn(xl * yl);
	double dxdz = 0;

	double dydx = dxdy;
	double dydy =  (n1*(n2-1)/(b*b)) * std::pow(std::abs(yl/b), n2 - 2) * std::pow(v, n1/n2 - 1) + (n2*(n1-n2)/(b*b)) * std::pow(std::abs(yl/b), 2*n2 - 2) * std::pow(v, n1/n2 - 2);
	double dydz = 0;
	
	double dzdx = 0;
	double dzdy = 0;
	double dzdz = (n1*(n1-1)/(c*c)) * std::pow(std::abs(zl/c), n1 - 2);

	Eigen::Matrix3d hess; 
	hess << dxdx, dxdy, dxdz, 
			dydx, dydy, dydz, 
			dzdx, dzdy, dzdz;
	return hess;
}


