#include "particle_generation.hpp"
#include "superellipsoid.hpp"
#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/src/Core/util/Constants.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <functional>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <random>

std::vector<Superellipsoid*>* generate_random_particles(std::vector<ParticleDistribution> particle_distributions, 
														std::vector<double> volume_fractions, 
														int n_particles) 
{
	if ( particle_distributions.size() != volume_fractions.size() ) {
		throw std::invalid_argument("[ERROR]: size of particle_distributions does not equal size of volume_fractions.");
	}

	// Checks that each particle distribution has unique cls value since having two
	// distribution for the same class of particle is not well defined.
	for (int ix = 0; ix < particle_distributions.size(); ix++) {
		int cls1 = particle_distributions.at(ix).cls;
		for (int jx = ix+1; jx < particle_distributions.size(); ix++) {
			int cls2 = particle_distributions.at(jx).cls;
			if ( cls1 == cls2 ) {
				throw std::invalid_argument("[ERROR]: each ParticleDistribution must have a unique cls value.");
			}
		}
	}

	// Compute sthe expected volumes for a particle from each article distribution
	// and computes the expected number of particles needed for each particle distribution
	// to on average get the target volume fraction 
	std::vector<double> exp_vols = expected_volumes(particle_distributions);
	std::vector<int> exp_n_particles = expected_particles_needed(exp_vols, volume_fractions);

	// Allocates vector on heap to store particles 
	std::vector<Superellipsoid*>* particles = new std::vector<Superellipsoid*>(n_particles);

	// Initialises random number generator.
	std::random_device rd;
	std::mt19937 mt(rd());
	
	// exp_vols, exp_n_particles, and particle_distributions are all index aligned
	// so we iterate over the particle distributions and generates the expected 
	// number of needed particles for each distribution and saves them to the heap
	int kx = 0; // Keeps track of index in particles vector 
	for (int ix = 0; particle_distributions.size(); ix++) {
		int cls = particle_distributions.at(ix).cls;
		Distribution d = particle_distributions.at(ix).volume_distribution;
		std::function<double(void)> sampler = get_sampler(d, mt);

		// Extracts parameter for reference particle
		Superellipsoid p = particle_distributions.at(ix).reference_particle;
		double scale_params[3] = {p.get_scale("a"), p.get_scale("b"), p.get_scale("c")};
		double shape_params[2] = {p.get_shape("n1"), p.get_shape("n2")};


		// TODO: enusre that sum(exp_n_particles) <= n_particles
		for (int jx = 0; jx < exp_n_particles.at(ix); jx++) {
			// Creates new superellipsoid on heap, scales it and 
			// put reference to it in return vector
			Superellipsoid* new_p = new Superellipsoid(cls, scale_params, shape_params);

		}
	}

	return particles;
}

Distribution parse_distribution(std::string distr_string) {
	int open_bracket_idx = distr_string.find("(");
	std::string distr_name = distr_string.substr(0, open_bracket_idx);
	std::string param_string = distr_string.substr(open_bracket_idx + 1, distr_string.size() - open_bracket_idx - 2);
	
	std::vector<double> param_vec;
	if (distr_name == "uniform") {
		int comma_idx = param_string.find(",");
		double lower_bound = std::stod(param_string.substr(0, comma_idx));
		double upper_bound = std::stod(param_string.substr(comma_idx + 1));
		
		param_vec.push_back(lower_bound);
		param_vec.push_back(lower_bound);

	} else if (distr_name == "normal") {
		int comma_idx = param_string.find(",");
		double mu = std::stod(param_string.substr(0, comma_idx));
		double sigma = std::stod(param_string.substr(comma_idx + 1));
		
		param_vec.push_back(mu);
		param_vec.push_back(sigma);

	} else if (distr_name == "log-normal") {
		int comma_idx = param_string.find(",");
		double mu = std::stod(param_string.substr(0, comma_idx));
		double sigma = std::stod(param_string.substr(comma_idx + 1));
			
		param_vec.push_back(mu);
		param_vec.push_back(sigma);

	} else {
		throw std::invalid_argument("[ERROR]: could not recognise distribution name. Available are: uniform, normal, log-normal.");
	}

	Distribution parsed_distr;
	parsed_distr.name = distr_name;
	parsed_distr.args = param_vec;
	return parsed_distr;
}

std::vector<int> expected_particles_needed(std::vector<double> expected_volumes, std::vector<double> target_volume_fractions) {
	// Allocates LHS and RHS of linear system
	Eigen::VectorXd b(expected_volumes.size());
	Eigen::MatrixXd A(expected_volumes.size(), expected_volumes.size());

	// Computes values for b and A entries
	for (int ix = 0; ix < A.rows(); ix++) {
		b(ix) = 0;
		for (int jx = 0; jx < A.cols(); jx++) {
			A(ix,jx) = (ix == jx) ? 
				expected_volumes.at(jx)*(1 - target_volume_fractions.at(jx)) : 
				expected_volumes.at(jx)*target_volume_fractions.at(jx);
		}
	}

	// Solves linear system using Householder rank-revealing QR decomp. with column pivoting. 
	Eigen::VectorXd x(expected_volumes.size());
	x = A.colPivHouseholderQr().solve(b);

	// Rounds and collects answer
	std::vector<int> n(expected_volumes.size());
	for (int ix = 0; ix < expected_volumes.size(); ix++) {
		n.at(ix) = (int) std::round(x(ix));
	}

	return n;
}

std::vector<double> expected_volumes(std::vector<ParticleDistribution> particle_distributions) {
	std::vector<double> ev(particle_distributions.size());
	for (int ix = 0; ix < particle_distributions.size(); ix++) {
		Distribution distr = particle_distributions.at(ix).volume_distribution;
		if (distr.name == "uniform") {
			double lb = distr.args.at(0); 
			double ub = distr.args.at(1);
			ev.at(ix) = (ub - lb)/2;
		} else if (distr.name == "normal") {
			ev.at(ix) = distr.args.at(0); // mu is at first index in vector
		} else if (distr.name == "lognormal") {
			// Our representation is on the form lognormal(mu, sigma)
			double mu = distr.args.at(0);
			double sigma = distr.args.at(1);
			ev.at(ix) = std::exp(mu + sigma*sigma/2);
		} else {
			throw std::invalid_argument("[ERROR]: could not recognise distribution name. Available are: uniform, normal, log-normal.");
		}
	}
	return ev;
}

std::function<double(void)> get_sampler(const Distribution& d, std::mt19937& mt) {
	// Checking so there exists exactly two arguments to d.
	// This obv. needs to be modified if more distributions are added
	if (d.args.size() != 2) throw std::invalid_argument("[ERROR]: the args vector to distribution must have length exactly size 2");

	if (d.name == "uniform") {
		std::uniform_real_distribution<double> sampler(d.args.at(0), d.args.at(1));
		return [sampler, &mt]()mutable{return sampler(mt);};
	} else if (d.name == "normal") {
		std::normal_distribution<double> sampler(d.args.at(0), d.args.at(1));
		return [sampler, &mt]()mutable{return sampler(mt);};
	} else if (d.name == "log-normal") {
		std::lognormal_distribution<double> sampler(d.args.at(0), d.args.at(1));
		return [sampler, &mt]()mutable{return sampler(mt);};
	} else {
		throw std::invalid_argument("[ERROR]: could not recognise distribution name. Available are: uniform, normal, log-normal.");
	}
}

Eigen::Quaternion<double> random_quaternion(std::mt19937& mt) {
	std::uniform_real_distribution<double> u(0.0, 1.0);
	double u1 = u(mt);
	double u2 = u(mt);
	double u3 = u(mt);

	Eigen::Quaternion<double> q(std::sqrt(1-u1) * std::sin(2 * M_PI * u2),
								std::sqrt(1-u1) * std::cos(2 * M_PI * u2),
								std::sqrt(u1) * std::sin(2 * M_PI * u3),
								std::sqrt(u1) * std::cos(2 * M_PI * u3));

	return q;
}



