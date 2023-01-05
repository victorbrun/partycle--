#include "particle_generation.hpp"
#include "superellipsoid.hpp" 
#include "coordinate_indexer.hpp" 
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

int main() {
	/*
	int n_particles = 100;
	std::vector<Superellipsoid*>* super_arr = new std::vector<Superellipsoid*>(n_particles);
	for (int ix = 0; ix < n_particles; ix++) {
		double scale_params[3] = {1,1,1};
		double shape_params[2] = {2,2};
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

	std::cout << "volume: " << super_arr->at(0)->volume() << std::endl;
	*/
	int n_particles = 100000;

	double ref_p1_scale[3] = {1,1,1};
	double ref_p2_scale[3] = {2,1,2};
	double ref_p3_scale[3] = {1,2,1};

	double ref_p1_shape[2] = {2,2};
	double ref_p2_shape[2] = {2,5};
	double ref_p3_shape[2] = {5,2};

	Superellipsoid ref_p1(1, ref_p1_scale, ref_p1_shape); 
	Superellipsoid ref_p2(2, ref_p2_scale, ref_p2_shape); 
	Superellipsoid ref_p3(3, ref_p3_scale, ref_p3_shape); 

	std::vector<double> ref_p1_vol_distr_args = {5, 7};
	std::vector<double> ref_p2_vol_distr_args = {10, 2};
	std::vector<double> ref_p3_vol_distr_args = {1, 0.25};

	std::vector<ParticleDistribution> pd = {
		{1, ref_p1, {"uniform", ref_p1_vol_distr_args}},
		{2, ref_p2, {"normal", ref_p2_vol_distr_args}},
		{3, ref_p3, {"log-normal", ref_p3_vol_distr_args}},
	};
	std::vector<double> target_volume_fractions = {0.2, 0.5, 0.3};

	// Generates random particles and saves them on the heap
	std::vector<Superellipsoid*>* particles = generate_random_particles(pd, 
																		target_volume_fractions,
																		n_particles);	

	return 0;
}
