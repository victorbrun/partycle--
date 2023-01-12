#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1	
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include "superellipsoid.hpp"

/***** Auxilirary functions *****/

// Sign function. Returns the sign of input and 0 if input is 0.
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/***** Private *****/

void Superellipsoid::set_scale(std::string param_name, double val) {
	if (param_name == "a") {
		this->scale[0] = val;	
	} else if (param_name == "b") {
		this->scale[1] = val;
	} else if (param_name == "c") {
		this->scale[2] = val;
	} else {
		throw std::invalid_argument("[ERROR]: invalid scale parameter name. Available are \"a\", \"b\", and \"c\"");
	}
}

/***** Public *****/

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

int Superellipsoid::get_class()  { return this->cls; }
std::vector<Superellipsoid*> get_contacts(Superellipsoid* p){return p->contacts; }
Eigen::Vector3d Superellipsoid::get_center()  { return this->c; }
Eigen::Quaternion<double> Superellipsoid::get_orientation()  { return this->q; }

void Superellipsoid::add_contact(Superellipsoid* p){this->contacts.push_back(p);}

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

double Superellipsoid::volume() {
	double a = this->get_scale("a");
	double b = this->get_scale("b");
	double c = this->get_scale("c");
	double n1 = this->get_shape("n1");
	double n2 = this->get_shape("n2");

	return 8 * a * b * c * (1/n1) * (1/n2) * std::beta(1/n1 + 1, 2/n1) * std::beta(1/n1, 1/n1);
}

double Superellipsoid::circumscribed_sphere_radius() {
	double a = this->get_scale("a");
	double b = this->get_scale("b");
	double c = this->get_scale("c");
	return std::sqrt( a*a + b*b + c*c );
}

double Superellipsoid::inscribed_sphere_radius() {
	double a = this->get_scale("a");
	double b = this->get_scale("b");
	double c = this->get_scale("c");
	return std::fmin( std::fmin(a, b), c );
}

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

void Superellipsoid::scale_to_volume(double vol) {
	double scale_factor = std::pow(vol/this->volume(), 1.0/3.0);
	
	double a = this->get_scale("a");
	double b = this->get_scale("b");
	double c = this->get_scale("c");

	this->set_scale("a", a*scale_factor);
	this->set_scale("b", b*scale_factor);
	this->set_scale("c", c*scale_factor);
}
Eigen::Vector2d Superellipsoid::constraints(Superellipsoid* p1, Superellipsoid* p2, Eigen::Vector4d& Z, double ax[3], double ay[3], double ex[2], double ey[2]) {
	
	std::vector<int> ind{0, 1, 2};
	Eigen::Vector3d X = Z(ind);

	// Particle parameters
	Eigen::Matrix3d Rx = p1->get_orientation().toRotationMatrix();
	Eigen::Matrix3d Ry = p2->get_orientation().toRotationMatrix();
	Eigen::Vector3d cx = p1->get_center();
	Eigen::Vector3d cy = p2->get_center();
	
	//Coordinates in local frames
    Eigen::Vector3d x = Rx.transpose() * (X-cx);
    Eigen::Vector3d y = Ry.transpose() * (X-cy);

	// Construct constraint
	Eigen::Vector2d H;
	H(0) = std::pow(std::pow(std::abs(x[0]/ax[0]), ex[1]) + std::pow(std::abs(x[1]/ax[1]), ex[1]), ex[0]/ex[1]) + std::pow(std::abs(x[2]/ax[2]), ex[0]) - 1;
    H(1) = std::pow(std::pow(std::abs(y[0]/ay[0]), ey[1]) + std::pow(std::abs(y[1]/ay[1]), ey[1]), ey[0]/ey[1]) + std::pow(std::abs(y[2]/ay[2]), ey[0]) - 1;

	return H;
}
Eigen::Matrix4d Superellipsoid::J_matrix(Superellipsoid* p1, Superellipsoid* p2, Eigen::Vector4d& Z, double ax[3], double ay[3], double ex[2], double ey[2]) {
	
	// Input :
	// 			Z : [X, mu]
	//			X : 3d global coordinates
	//			mu: double 
	std::vector<int> ind{0, 1, 2};
	Eigen::Vector3d X = Z(ind);
	double mu = Z(3);

	// Particle parameters
	Eigen::Matrix3d Rx = p1->get_orientation().toRotationMatrix();
	Eigen::Matrix3d Ry = p2->get_orientation().toRotationMatrix();
	Eigen::Vector3d cx = p1->get_center();
	Eigen::Vector3d cy = p2->get_center();

	//Coordinates in local frames
    Eigen::Vector3d x = Rx.transpose() * (X-cx);
    Eigen::Vector3d y = Ry.transpose() * (X-cy);

	// Temporary variables
    double vx = std::pow(std::abs(x[0]/ax[0]), ex[1]) + std::pow(std::abs(x[1]/ax[1]), ex[1]);
    double vy = std::pow(std::abs(y[0]/ay[0]), ey[1]) + std::pow(std::abs(y[1]/ay[1]), ey[1]);
    double k1 = ex[0]/ex[1]-1;
    double k2 = ey[0]/ey[1]-1;

    // Local first derivatives
    double f1x = ex[0]/ax[0] * std::pow(std::abs(x[0]/ax[0]), (ex[1]-1)) * std::pow(vx, k1) * sgn(x[0]);
    double f1y = ex[0]/ax[1] * std::pow(std::abs(x[1]/ax[1]), (ex[1]-1)) * std::pow(vx, k1) * sgn(x[1]);        
    double f1z = ex[0]/ax[2] * std::pow(std::abs(x[2]/ax[2]), (ex[0]-1)) * sgn(x[2]);

	double f2x = ey[0]/ay[0] * std::pow(std::abs(y[0]/ay[0]), (ey[1]-1)) * std::pow(vy, k2) * sgn(y[0]);
    double f2y = ey[0]/ay[1] * std::pow(std::abs(y[1]/ay[1]), (ey[1]-1)) * std::pow(vy, k2) * sgn(y[1]);   
    double f2z = ey[0]/ay[2] * std::pow(std::abs(y[2]/ay[2]), (ey[0]-1)) * sgn(y[2]);

	// Local gradients
    Eigen::Vector3d f1g {f1x, f1y, f1z}; 
    Eigen::Vector3d f2g {f2x, f2y, f2z};  

    // Global gradients
    Eigen::Vector3d F1g = Rx * f1g;
    Eigen::Vector3d F2g = Ry * f2g;

    // Constraint hessians
    Eigen::Matrix3d h1 = p1->inside_outside_hess(x);
    Eigen::Matrix3d h2 = p2->inside_outside_hess(x);

    Eigen::Matrix3d H1 = Rx * h1 * Rx.transpose();
    Eigen::Matrix3d H2 = Ry * h2 * Ry.transpose();


    // Construct J
	Eigen::Matrix4d J;

    J(ind, ind) = H1 + std::pow(mu, 2) * H2;
    J(3  , ind) = F1g-F2g;
    J(ind,   3) = 2*mu*F2g;

    return J;
}
Eigen::Vector4d Superellipsoid::phi(Superellipsoid* p1, Superellipsoid* p2, Eigen::Vector4d &Z, double ax[3], double ay[3], double ex[2], double ey[2]) {

	// Input :	Z : [X, mu]
	//			X : 3d global coordinates
	//			mu: double 
	std::vector<int> ind{0, 1, 2};
	Eigen::Vector3d X = Z(ind);
	double mu = Z(3);

	// Particle parameters
	Eigen::Matrix3d Rx = p1->get_orientation().toRotationMatrix();
	Eigen::Matrix3d Ry = p2->get_orientation().toRotationMatrix();
	Eigen::Vector3d cx = p1->get_center();
	Eigen::Vector3d cy = p2->get_center();

	//Coordinates in local frames
    Eigen::Vector3d x = Rx.transpose() * (X-cx);
    Eigen::Vector3d y = Ry.transpose() * (X-cy);

	// Temporary variables
    double vx = std::pow(std::abs(x[0]/ax[0]), ex[1]) + std::pow(std::abs(x[1]/ax[1]), ex[1]);
    double vy = std::pow(std::abs(y[0]/ay[0]), ey[1]) + std::pow(std::abs(y[1]/ay[1]), ey[1]);
    double k1 = ex[0]/ex[1]-1;
    double k2 = ey[0]/ey[1]-1;
 	
    // Local first derivatives
    double f1x = ex[0]/ax[0] * std::pow(std::abs(x[0]/ax[0]), (ex[1]-1)) * std::pow(vx, k1) * sgn(x[0]);
    double f1y = ex[0]/ax[1] * std::pow(std::abs(x[1]/ax[1]), (ex[1]-1)) * std::pow(vx, k1) * sgn(x[1]);       
    double f1z = ex[0]/ax[2] * std::pow(std::abs(x[2]/ax[2]), (ex[0]-1)) * sgn(x[2]);

	double f2x = ey[0]/ay[0] * std::pow(std::abs(y[0]/ay[0]), (ey[1]-1)) * std::pow(vy, k2) * sgn(y[0]);
    double f2y = ey[0]/ay[1] * std::pow(std::abs(y[1]/ay[1]), (ey[1]-1)) * std::pow(vy, k2) * sgn(y[1]);      
    double f2z = ey[0]/ay[2] * std::pow(std::abs(y[2]/ay[2]), (ey[0]-1)) * sgn(y[2]);

	// Local gradients
    Eigen::Vector3d f1g {f1x, f1y, f1z}; 
    Eigen::Vector3d f2g {f2x, f2y, f2z}; 

    // Global gradients
    Eigen::Vector3d F1g = Rx * f1g;
    Eigen::Vector3d F2g = Ry * f2g;
	
    double F1 = std::pow(std::pow(std::abs(x[0]/ax[0]), ex[1]) + std::pow(std::abs(x[1]/ax[1]), ex[1]), ex[0]/ex[1]) + std::pow(std::abs(x[2]/ax[2]), ex[0]) - 1;
    double F2 = std::pow(std::pow(std::abs(y[0]/ay[0]), ey[1]) + std::pow(std::abs(y[1]/ay[1]), ey[1]), ey[0]/ey[1]) + std::pow(std::abs(y[2]/ay[2]), ey[0]) - 1;
    
	// Construct Phi
	Eigen::Vector4d phi;
    phi(ind) = F1g + std::pow(mu, 2) * F2g;
    phi(3) = F1 - F2;

	return phi;
}

bool Superellipsoid::distance(Superellipsoid* p1, Superellipsoid* p2){
	// Parameters
	int N = 4;
	int K = 8;
	double eps = 1e-1;

	// Particle parameters

	double ax[3] = {p1->scale[0], p1->scale[1], p1->scale[2]};
	double ay[3] = {p2->scale[0], p2->scale[1], p2->scale[2]};
	double ex[2] = {p1->shape[0], p1->shape[1]};
	double ey[2] = {p2->shape[0], p2->shape[1]};
	Eigen::Vector3d cx = p1->get_center();
	Eigen::Vector3d cy = p2->get_center();
	// Smallest scale parameters
	Eigen::Vector3d tx = {ax[0], ax[1], ax[2]};
	Eigen::Vector3d ty = {ay[0], ay[1], ay[2]};
	double rx = tx.minCoeff(); 
	double ry = ty.minCoeff();

	// Initial and step of shape and scale parameters
	double dax[3] =  {(ax[0]-rx)/N, (ax[2]-rx)/N, (ax[2]-rx)/N};
	double day[3] =  {(ay[0]-ry)/N, (ay[2]-ry)/N, (ay[2]-ry)/N};
	double dex[2] =  {(ex[0]-2)/N, (ex[1]-2)/N};
	double dey[2] =  {(ey[0]-2)/N, (ey[1]-2)/N};
	double ax0[3] =  {rx, rx, rx};
	double ay0[3] =  {ry, ry, ry};
	double ex0[2] =  {2, 2};
	double ey0[2] =  {2, 2};

	// Initialize coordinates to x0, 1
	double pNi; 
	double pNo;
	Eigen::Vector4d Z = {(cx[0]+cy[0])/2, (cx[1]+cy[1])/2, (cx[2]+cy[2])/2, 1};
	bool b;
	for (int i = 0; i < N; i++) {
		// Increment shape and scale parameters
		double exi[2] = {ex0[0] + i * dex[0], ex0[1] + i * dex[1]};
		double eyi[2] = {ey0[0] + i * dey[0], ey0[1] + i * dey[1]};
        double axi[3] = {ax0[0] + i * dax[0], ax0[1] + i * dax[1], ax0[2] + i * dax[2]};
		double ayi[3] = {ay0[0] + i * day[0], ay0[1] + i * day[1], ay0[2] + i * day[2]};
		
		for (int j = 0;  j< K; j++) {
			// For select shape and scale parameters, iterate solution

			// Collect phi and J
			Eigen::Vector4d phi = Superellipsoid::phi(p1, p2, Z, axi, ayi, exi, eyi);
            Eigen::Matrix4d J   = Superellipsoid::J_matrix(p1, p2, Z, axi, ayi, exi, eyi);

			// Solve equation system : 
			// TODO: Test different solvers for speed
			Eigen::Vector4d dZ = J.colPivHouseholderQr().solve(-phi);
			Eigen::Vector4d Zp;
			pNo = phi.norm();
			double alpha = 1;
			double d;
			while (true){
				Zp = Z+alpha*dZ;
				pNi = Superellipsoid::phi(p1, p2, Zp, axi, ayi, exi, eyi).norm();
				if (pNi<pNo) {
					// If we're closer to a solution; accept step size
					d = alpha*dZ.norm();
					Z = Z+alpha*dZ;
					break;

				} else {
					alpha = alpha/2;

				}
			if(d<eps) {
				break;
			}
			}
		}
	
	}
    b = false;
	Eigen::Vector2d h = Superellipsoid::constraints(p1, p2, Z, ax, ay, ex, ey);
    if (h(1)<0 && h(1)<0) {
		// Collision
        b = true;
	}
        
    return b;
}

