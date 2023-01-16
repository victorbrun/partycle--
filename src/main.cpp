#include "particle_generation.hpp"
#include "program_options.hpp"
#include "superellipsoid.hpp"
#include "utils.hpp"
#include "domain.hpp"

#include <algorithm>
#include <chrono>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#define CONTACT_TOL 1e-6

int main(int argc, char* argv[]) {
	// Parses input arguments and stores them in file-scope varaibel	
	program_options::parse(argc, argv);

	// Extract the contact tolerance from parsed input. Note that
	// there is no catch for when input tolerance cannot be converted to double.
	double contact_tol;
	if (program_options::get_option("--contact-tolerance") != "") {
		contact_tol = std::stod(program_options::get_option("--contact-tolerance"));
	} else if (program_options::get_option("-ct") != "") {
		contact_tol = std::stod(program_options::get_option("-ct"));
	} else {
		contact_tol = CONTACT_TOL;
	}

	// Extracts domain from parsed input 
	std::string domain_string;
	if (program_options::get_option("--domain") != "") {
		domain_string = program_options::get_option("--domain");
	} else if (program_options::get_option("-d") != "") {
		domain_string = program_options::get_option("-d");
	} else {
		std::cerr << "partycle--: need domain to be defined" << "\n";
		std::cerr << "usage: partycle-- [-d|--domain] [ax,bx]x[ay,by]x[az,bz] [-cf|--component-file] <input_file>...\n";
		return EXIT_FAILURE;
	}

	// Tries to parse the domain string into bounds for the domain 
	std::vector<double> domain_bounds;
	try {
		domain_bounds = parse_domain(domain_string);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << "\n";
		return EXIT_FAILURE;
	}
	
	// Extracts the file containing the information about each component in the mixture
	std::string components_file_name;
	if (program_options::get_option("--component-file") != "") {
		components_file_name = program_options::get_option("--component-file");
	} else if (program_options::get_option("-cf") != "") {
		components_file_name = program_options::get_option("-cf");
	} else {
		std::cerr << "partycle--: CSV file containing mixture components must be provided" << "\n";
		std::cerr << "usage: partycle-- [-d|--domain] [ax,bx]x[ay,by]x[az,bz] [-cf|--component-file] <input_file>...\n";
		return EXIT_FAILURE;
	}

	// Tries to parse the contents of the component-file
	std::vector<Component> components;
	try {
		components = parse_components(components_file_name);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << "\n";
		return EXIT_FAILURE;
	} catch (const std::invalid_argument& err) {
		std::cerr << err.what() << "\n";
		return EXIT_FAILURE;
	}

	// Defines domain bounds 
	double x_range[2] = {domain_bounds.at(0), domain_bounds.at(1)};
	double y_range[2] = {domain_bounds.at(2), domain_bounds.at(3)};
	double z_range[2] = {domain_bounds.at(4), domain_bounds.at(5)};
	double domain_vol = (x_range[1] - x_range[0]) * (y_range[1] - y_range[0]) * (z_range[1] - z_range[0]);

	// Generates particles acocrding to the specification given in component-file
	std::vector<Superellipsoid*>* particles = generate_random_particles(components, domain_vol);
	std::cout << "[INFO]: " << particles->size() << " particles generated" << std::endl;

	// Shuffles the particles so that by iterating over the vector 
	// will be the same as randomly selecting a component and generating 
	// a particle according to it. 
	// TODO: check that this argument of shuffling pre-generated particles is equiv. to that of randomly sampling
	// component and then generating particle.
	std::random_device rd;
	std::mt19937 mt(rd());
	std::shuffle(particles->begin(), particles->end(), mt);

	// Initialises domain. This is done here and not where its bounds are defined since we want to 
	// specify the number of particles we are going to add to it. This will make the domain pre-allocate
	// memory.
	Domain domain = Domain(x_range, y_range, z_range, contact_tol, particles->size());
	
	// Uses first four particles to initiate advancing front
	Superellipsoid* init_p[4] = {particles->at(0), particles->at(1), particles->at(2), particles->at(3)};
	domain.initialise_outward_advancing_front(init_p);
	
	// Uses the rest of the particles to incerement the advanding front
	for (size_t ix = 4; ix < particles->size(); ix++) {
		//std::cout << "[INFO]: placing particle number: " << ix+1 << std::endl;
		//auto start = std::chrono::high_resolution_clock::now();
		domain.increment_advancing_front(particles->at(ix));
		//auto stop = std::chrono::high_resolution_clock::now();
		//std::cout << "total time to place particle " << ix+1 << ": " << (double(std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()))/1000000.0 << "s" << std::endl;
	}
	domain.write_csv("domain.csv");

	// TODO: fill the domain using advancing front!!
	// TODO: Compute contact statistics and output it in some reasonable way
	// TODO: DONE!

	std::cout << "[INFO]: program finished, cleaning up.. " << std::endl;
	std::cout << "[INFO]: destroying generated particles" << std::endl;
	delete particles;	

	return EXIT_SUCCESS;
}
