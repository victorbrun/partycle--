#ifndef SUPERELLIPSOID_H
#define SUPERELLIPSOID_H

#include <string>
#include <eigen3/Eigen/Dense>

class Superellipsoid {
	private:
		int 						cls;
		Eigen::Vector3d 			c;
		double 						scale[3];
		double 						shape[2];
		Eigen::Quaternion<double> 	q;

		Eigen::Vector3d to_local_coords(Eigen::Vector3d x);

	public:
		Superellipsoid(int cls, double scale_params[3], double shape_params[2]);

		// Getters
		int get_class();
		Eigen::Vector3d get_center();
		double get_scale(std::string param_name);
		double get_shape(std::string param_name);
		Eigen::Quaternion<double> get_orientation();

		// Setters
		void set_center(Eigen::Vector3d c);
		void set_orientation(Eigen::Quaternion<double> q);
	
		double implicit_surface(Eigen::Vector3d x);
		double inside_outside(Eigen::Vector3d x);
		Eigen::Vector3d inside_outside_grad(Eigen::Vector3d x);
		Eigen::Matrix3d inside_outside_hess(Eigen::Vector3d x);

		double parametric_surface(double eta, double omega);
		double enclosing_sphere_radius();
		void draw(int resolution);
		static double distance(Superellipsoid p1, Superellipsoid p2);

};

#endif
