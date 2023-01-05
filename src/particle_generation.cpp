#include "particle_generation.hpp"
#include "superellipsoid.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstddef>
#include <eigen3/Eigen/src/Core/util/Constants.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <functional>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <random>

std::vector<Superellipsoid*>* generate_random_particles(const std::vector<ParticleDistribution>& particle_distributions, 
														const std::vector<double>& volume_fractions,
														double domain_volume) 
{
	if ( particle_distributions.size() != volume_fractions.size() ) {
		throw std::invalid_argument("[ERROR]: size of particle_distributions does not equal size of volume_fractions.");
	}

	// Checks that each particle distribution has unique cls value since having two
	// distribution for the same class of particle is not well defined.
	for (int ix = 0; ix < particle_distributions.size(); ix++) {
		int cls1 = particle_distributions.at(ix).cls;
		for (int jx = ix+1; jx < particle_distributions.size(); jx++) {
			int cls2 = particle_distributions.at(jx).cls;
			if ( cls1 == cls2 ) {
				throw std::invalid_argument("[ERROR]: each ParticleDistribution must have a unique cls value.");
			}
		}
	}
	// Initialises random number generator.
	std::random_device rd;
	std::mt19937 mt(rd());

	// Computes the expected number of particles needed of each class to achieve `volume_fractions`
	std::vector<int> exp_n_particles = expected_particles_needed(particle_distributions, volume_fractions, domain_volume);
	int total_n_particles = std::accumulate(exp_n_particles.begin(), exp_n_particles.end(), 0);

	// Allocates vector on heap to store particles 
	std::vector<Superellipsoid*>* particles = new std::vector<Superellipsoid*>(total_n_particles);
	
	size_t kx = 0; // Keeps track of next index in particles
	for (size_t ix = 0; ix < exp_n_particles.size(); ix++) {
		ParticleDistribution pd = particle_distributions.at(ix);
		std::function<double(void)> volume_sampler = get_sampler(pd.volume_distribution, mt);

		// Extracts parameters for reference particle 
		Superellipsoid p = pd.reference_particle;
		int cls = pd.cls;
		double scale_params[3] = {p.get_scale("a"), p.get_scale("b"), p.get_scale("c")};
		double shape_params[2] = {p.get_shape("n1"), p.get_shape("n2")};
		
		// Genereates new particles with random volume and rotation
		for (size_t jx = 0; jx < exp_n_particles.at(ix); jx++) {
			Superellipsoid* new_p = new Superellipsoid(cls, scale_params, shape_params);
			new_p->set_orientation(random_quaternion(mt));
			new_p->scale_to_volume(volume_sampler());
			particles->at(kx) = new_p;
			kx++;
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

std::vector<int> expected_particles_needed(const std::vector<ParticleDistribution>& particle_distributions, 
										   const std::vector<double>& target_volume_fractions,
										   double domain_volume) {
	std::vector<int> n(particle_distributions.size());
	for (size_t ix = 0; ix < n.size(); ix++) {
		Distribution d = particle_distributions.at(ix).volume_distribution;
		double eta = target_volume_fractions.at(ix);

		n.at(ix) = eta*domain_volume/mean(d);
	}

	return n;
}

std::vector<double> particle_class_selection_prob(const std::vector<ParticleDistribution>& particle_distributions,
												  const std::vector<double>& target_volume_fractions) {
	std::vector<double> p(particle_distributions.size());
	
	double mean_sum = 0;
	for (size_t ix = 0; ix < p.size(); ix++) {
		mean_sum += mean(particle_distributions.at(ix).volume_distribution);
	}

	for (size_t ix = 0; ix < p.size(); ix++) {
		double mu = mean(particle_distributions.at(ix).volume_distribution);
		double eta = target_volume_fractions.at(ix);
		p.at(ix) = ( p.size() * mu * eta )/mean_sum; 
	}

	return p;
}

double mean(const Distribution& d) {
	if (d.name == "uniform") {
		double lb = d.args.at(0); 
		double ub = d.args.at(1);
		return (ub - lb)/2;
	} else if (d.name == "normal") {
		return d.args.at(0); // mu is at first index in vector
	} else if (d.name == "log-normal") {
		// Our representation is on the form lognormal(mu, sigma)
		double mu = d.args.at(0);
		double sigma = d.args.at(1);
		return std::exp(mu + sigma*sigma/2);
	} else {
		throw std::invalid_argument("[ERROR]: could not recognise distribution name. Available are: uniform, normal, log-normal.");
	}
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



