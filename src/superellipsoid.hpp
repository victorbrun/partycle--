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

		// Changes the scale parameter of Superellipsoid to `val`
		void set_scale(std::string param_name, double val);

		static Eigen::Matrix4d J_matrix(const Superellipsoid& p1, const Superellipsoid& p2, const Eigen::Vector3d& x);

		static Eigen::Vector3d phi_vec(const Superellipsoid& p1, const Superellipsoid& p2, const Eigen::Vector3d& x);

	public:
		Superellipsoid(int cls, double scale_params[3], double shape_params[2]);

		// Returns class of Superellipsoid.
		int get_class();

		// Returns center of Superellipsoid
		Eigen::Vector3d get_center();

		/**
		 * Returns a scale parameter for Superellipsoid.
		 *
		 * @param `param_name`: name of scale parameter. Available are: a, b, c.
		 * @return scale parameter value.
		 */
		double get_scale(std::string param_name);
		
		/**
		 * Returns a shape parameter for Superellipsoid.
		 *
		 * @param `param_name`: name of shape parameter. Available are: n1, n2.
		 * @return shape parameter value.
		 */
		double get_shape(std::string param_name);
		
		//Returns the orientation of the Superellipsoid as a Quaternion.
		Eigen::Quaternion<double> get_orientation();

		// Sets the center of the Superellipsoid
		void set_center(Eigen::Vector3d c);

		// Sets the orientation of the Superellipsoid
		void set_orientation(Eigen::Quaternion<double> q);
	
		// Computes the volume for the current Superellipsoid
		double volume();
	
		/*
		 * Computes the value of left hand side of below implicit function, which defines the Superellipsoid surface:
		 * 		( |x/a|^n2 + |y/b|^n2 )^(n1/n2) + |z/c|^n1 = 1
		 *
		 * @param `x`: vector in R^3. `x = (x, y, z)^T`.
		 * @return value of above LHS evaluated at `x`.
		 */
		double implicit_surface(Eigen::Vector3d x);

		/**
		 * Implicit surface funciton -1. Can be used to indicate if point is inside, on or outside Superellipsoid.
		 *
		 * @param `x`: vector in R^3. `x = (x, y, z)^T`.
		 * @return value of implicit surface function evaluate ad `x` -1. 
		 */
		double inside_outside(Eigen::Vector3d x);

		// Gradient vector of inside_outside function evaluated at `x = (x, y, z)^T`.
		Eigen::Vector3d inside_outside_grad(Eigen::Vector3d x);

		// Hessian matrix of inside_outside function evalueate at `x = (x, y, z)^T`.
		Eigen::Matrix3d inside_outside_hess(Eigen::Vector3d x);

		/*
		 * Computes the parametric radius vector of the Superellipsoid as ´eta´, ´omega´.
		 *
		 * @param `eta`: space angle, defined on [-pi/2, pi/2].
		 * @param `omega`: xy-place angle, defined on [-pi, pi).
		 * @return radius vector evaluated at `eta`, `omega`.
		 */
		double parametric_surface(double eta, double omega);
		
		// Returns radius of smallest sphere enclosing whole Superellipsoid.
		double circumscribed_sphere_radius();

		// Returns radius of largest sphere enclosed by Superellipsoid
		double inscribed_sphere_radius();
		
		/**
		 * Returns bool indicating if computation succeeded and minimum 
		 * distance between `p1` and `p2`.
		 */		
		static bool cc(Superellipsoid* p1, Superellipsoid* p2);

		// Re-scales particle to have the volume specified by `vol`.
		void scale_to_volume(double vol);
};

#endif
