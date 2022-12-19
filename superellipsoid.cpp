#include <boost/math/quaternion.hpp>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <stdexcept>
#include <string>
#include "superellipsoid.h"

Superellipsoid::Superellipsoid(int cls, float scale_params[3], float shape_params[2]) {
	float n1 = shape_params[0];
	float n2 = shape_params[1];
	if (n1 < 2 || n2 < 2) {
		throw std::invalid_argument("ERROR: n1 and n2 must be larger or equal to 2 to ensure that superellipsoid surface is two times differentiable");
	}
	this->shape[0] = n1;
	this->shape[1] = n2;

	float a = scale_params[0];
	float b = scale_params[1];
	float c = scale_params[2];
	this->scale[0] = a;
	this->scale[1] = b;
	this->scale[2] = c;
		
	this->cls = cls;
}

int Superellipsoid::get_class() { return this->cls; }
Eigen::Vector3d Superellipsoid::get_center() { return this->c; }
Eigen::Quaternion<double> Superellipsoid::get_orientation() { return this->q; }

float Superellipsoid::get_scale(std::string param_name) {
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

float Superellipsoid::get_shape(std::string param_name) {
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

float Superellipsoid::implicit_surface(Eigen::Vector3d x) {
	Eigen::Matrix3d R = this->get_orientation().toRotationMatrix();
	Eigen::Vector3d x_local =  R * ( x - this->get_center() );
	float xl = x_local[0];
	float yl = x_local[1];
	float zl = x_local[2];

	float a 	= this->get_scale("a"); 
	float b 	= this->get_scale("b"); 
	float c 	= this->get_scale("c");
	float n1 	= this->get_shape("n1");
	float n2 	= this->get_shape("n2");

	return std::pow( std::pow(std::abs(xl/a), n2) + std::pow(std::abs(yl/b), n2), n1/n2 ) + std::pow(std::abs(zl/c), n1);

	return 0;
}

float Superellipsoid::inside_outside(Eigen::Vector3d x) {
	return this->implicit_surface(x) - 1;
}



