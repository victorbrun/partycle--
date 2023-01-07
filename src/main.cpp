#include "particle_generation.hpp"
#include "superellipsoid.hpp" 
#include "coordinate_indexer.hpp" 
#include "domain.hpp"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>
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

	// Domain is to large for total volume to be stored in double so it overflows
	double x_range[2] = {0,100};
	double y_range[2] = {0,100};
	double z_range[2] = {0,100};
	double domain_vol = (x_range[1]-x_range[0])*(y_range[1]-y_range[0])*(z_range[1]-z_range[0]); 
	std::vector<int> exp_n_particles = expected_particles_needed(pd, target_volume_fractions, domain_vol);
	int n_particles = std::accumulate(exp_n_particles.begin(), exp_n_particles.end(), 0);
	Domain d = Domain(x_range, y_range, z_range, n_particles);
	
	// Generates random particles and saves them on the heap
	size_t runs = 50;
	Eigen::ArrayXXd samples = Eigen::ArrayXXd::Zero(3, runs);
	std::cout << "Generating " << n_particles << " random particles of " << pd.size() << " different classes " << runs << " times." <<std::endl;
	for (size_t jx = 0; jx < runs; jx++) {		
		std::vector<Superellipsoid*>* particles = generate_random_particles(pd, target_volume_fractions, domain_vol);
		for (size_t ix = 0; ix < particles->size(); ix++) {
			Superellipsoid* p = particles->at(ix);
			samples(p->get_class()-1, jx) += p->volume();
			//std::cout << "samples(" << p->get_class()-1 << "," << jx << ") = " << samples(p->get_class()-1, jx) << std::endl;
		}
		delete particles;
	}

	// Computes average volume over above runns
	long double average_vol_arr[3] = {0, 0, 0};
	for (size_t ix = 0; ix < runs; ix++) {
		average_vol_arr[0] += samples(0, ix)/runs;
		average_vol_arr[1] += samples(1, ix)/runs;
		average_vol_arr[2] += samples(2, ix)/runs;

		//std::cout << average_vol_arr[0] << " " << average_vol_arr[1] << " " << average_vol_arr[2] << std::endl;
	}
	

	
	long double vol_sum = average_vol_arr[0] + average_vol_arr[1] + average_vol_arr[2];
	std::cout << " ---------- SEQUENTIALLY GENERATED PARTICLES -------------" << std::endl;
	std::cout << "[INFO]: Target volume fractions: class 1: " << target_volume_fractions.at(0) << " class 2: " << target_volume_fractions.at(1) << " class 3: " << target_volume_fractions.at(2) << std::endl;  
	std::cout << "[INFO]: Actual volume fractions: class 1: " << average_vol_arr[0]/vol_sum << " class 2: " << average_vol_arr[1]/vol_sum << " class 3: " << average_vol_arr[2]/vol_sum << std::endl;  
	std::cout << "\n";	

	return 0;
}
