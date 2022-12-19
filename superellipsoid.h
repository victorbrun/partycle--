#ifndef SUPERELLIPSOID_H
#define SUPERELLIPSOID_H

#include <string>
#include <eigen3/Eigen/Dense>

class Superellipsoid {
	private:
		int 						cls;
		Eigen::Vector3d 			c;
		float 						scale[3];
		float 						shape[2];
		Eigen::Quaternion<double> 	q;


	public:
		Superellipsoid(int cls, float scale_params[3], float shape_params[2]);

		// Getters
		int get_class();
		Eigen::Vector3d get_center();
		float get_scale(std::string param_name);
		float get_shape(std::string param_name);
		Eigen::Quaternion<double> get_orientation();

		// Setters
		void set_center(Eigen::Vector3d c);
		void set_orientation(Eigen::Quaternion<double> q);
	
		float implicit_surface(Eigen::Vector3d x);
		float inside_outside(Eigen::Vector3d x);
		Eigen::Vector3d inside_outside_grad(Eigen::Vector3d x);
		Eigen::Matrix3d inside_outside_hess(Eigen::Vector3d x);

		float parametric_surface(float eta, float omega);
		float enclosing_sphere_radius();
		void draw(int resolution);
		static float distance(Superellipsoid p1, Superellipsoid p2);
};

#endif
