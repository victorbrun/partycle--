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

std::vector<Superellipsoid*>* generate_random_particles(const std::vector<Component>& components, const double domain_volume) 
{
	// Checks that each component has unique component_id value since having two volume
	// distributions for the same component is not well defined (yet..). Also checks that
	// the target volume fractions sum to one.
	double frac_sum = 0;
	for (size_t ix = 0; ix < components.size(); ix++) {
		int cls1 = components.at(ix).component_id;
		for (size_t jx = ix+1; jx < components.size(); jx++) {
			int cls2 = components.at(jx).component_id;
			if ( cls1 == cls2 ) {
				throw std::invalid_argument("[ERROR]: each ParticleDistribution must have a unique cls value.");
			}
		}
		frac_sum += components.at(ix).target_volume_fraction;
	}
	// Initialises random number generator.
	std::random_device rd;
	std::mt19937 mt(rd());

	// Computes the expected number of particles needed of each class to achieve `volume_fractions`
	std::vector<int> exp_n_particles = expected_particles_needed(components, domain_volume);
	int total_n_particles = std::accumulate(exp_n_particles.begin(), exp_n_particles.end(), 0);
	
	std::cout << "[INFO]: generating " << total_n_particles << " particles distributed over the components: ";  
	for (size_t ix = 0; ix < exp_n_particles.size(); ix++) {
		std::cout << components.at(ix).component_id << ": " << exp_n_particles.at(ix) << ", ";
	}
	std::cout << "\n";

	// Allocates vector on heap to store particles 
	std::vector<Superellipsoid*>* particles = new std::vector<Superellipsoid*>(total_n_particles);
	
	size_t kx = 0; // Keeps track of next index in particles
	for (size_t ix = 0; ix < exp_n_particles.size(); ix++) {
		Component cmp = components.at(ix);
		std::function<double(void)> volume_sampler = get_sampler(cmp.volume_distribution, mt);

		// Extracts parameters for reference particle 
		Superellipsoid p = cmp.reference_particle;
		int component_id = cmp.component_id;
		double scale_params[3] = {p.get_scale("a"), p.get_scale("b"), p.get_scale("c")};
		double shape_params[2] = {p.get_shape("n1"), p.get_shape("n2")};
		
		// Genereates new particles with random volume and rotation
		for (size_t jx = 0; jx < (size_t)exp_n_particles.at(ix); jx++) {
			Superellipsoid* new_p = new Superellipsoid(component_id, scale_params, shape_params);
			new_p->set_orientation(Eigen::Quaternion<double>::UnitRandom());
			new_p->scale_to_volume(volume_sampler());
			particles->at(kx) = new_p;
			kx++;
		}
	}

	return particles;
}

Distribution parse_distribution(const std::string& distr_string) {
	int open_bracket_idx = distr_string.find("(");
	std::string distr_name = distr_string.substr(0, open_bracket_idx);
	std::string param_string = distr_string.substr(open_bracket_idx + 1, distr_string.size() - open_bracket_idx - 2);
	
	std::vector<double> param_vec;
	if (distr_name == "uniform") {
		int comma_idx = param_string.find(",");
		double lower_bound = std::stod(param_string.substr(0, comma_idx));
		double upper_bound = std::stod(param_string.substr(comma_idx + 1));

		// Making sure distribution arguments make sense
		if (lower_bound > upper_bound) {
			throw std::invalid_argument("[ERROR]: first argument must be smaller than second argument in uniform distribution.");
		}

		param_vec.push_back(lower_bound);
		param_vec.push_back(upper_bound);

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

	} else if (distr_name == "weibull") {
		// Our representation is on the form weibull(k, lambda)
		int comma_idx = param_string.find(",");
		double k = std::stod(param_string.substr(0, comma_idx)); // shape parameter 
		double lambda = std::stod(param_string.substr(comma_idx + 1)); // scale parameter
			
		param_vec.push_back(k);
		param_vec.push_back(lambda);
	} else {
		throw std::invalid_argument("[ERROR]: could not recognise distribution name. Available are: uniform, normal, log-normal.");
	}

	Distribution parsed_distr;
	parsed_distr.name = distr_name;
	parsed_distr.args = param_vec;
	return parsed_distr;
}

std::vector<int> expected_particles_needed(const std::vector<Component>& components, const double domain_volume) {
	std::vector<int> n(components.size());
	for (size_t ix = 0; ix < n.size(); ix++) {
		Distribution d = components.at(ix).volume_distribution;
		double eta = components.at(ix).target_volume_fraction;

		n.at(ix) = eta*domain_volume/mean(d);
	}

	return n;
}

std::vector<double> particle_class_selection_prob(const std::vector<Component>& components,
												  const double domain_volume,
												  const size_t n_particles) {
	std::vector<double> p(components.size());
	
	for (size_t ix = 0; ix < p.size(); ix++) {
		double f = components.at(ix).target_volume_fraction;
		Distribution vol_d = components.at(ix).volume_distribution;

		p.at(ix) = domain_volume * f /(n_particles * mean(vol_d));
	}

	return p;
}

double mean(const Distribution& d) {
	if (d.name == "uniform") {
		double lb = d.args.at(0); 
		double ub = d.args.at(1);
		return (lb + ub)/2;
	} else if (d.name == "normal") {
		return d.args.at(0); // mu is at first index in vector
	} else if (d.name == "log-normal") {
		// Our representation is on the form lognormal(mu, sigma)
		double mu = d.args.at(0);
		double sigma = d.args.at(1);
		return std::exp(mu + sigma*sigma/2);
	} else if (d.name == "weibull") {
		// Our representation is on the form weibull(k, lambda)
		double k = d.args.at(0);
		double lambda = d.args.at(1);
		return lambda * std::tgamma(1 + 1/k);
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
	} else if (d.name == "weibull") {
		// Our representation is on the form weibull(k, lambda)
		std::weibull_distribution<double> sampler(d.args.at(0), d.args.at(1));
		return [sampler, &mt]()mutable{return sampler(mt);};
	} else {
		throw std::invalid_argument("[ERROR]: could not recognise distribution name. Available are: uniform, normal, log-normal.");
	}
}