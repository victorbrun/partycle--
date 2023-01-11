#include "particle_generation.hpp"
#include "program_options.hpp"
#include "superellipsoid.hpp"
#include "utils.hpp"
#include "domain.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#define CONTACT_TOL 1e-3

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

	// Creates domain in which to place particles 
	double x_range[2] = {domain_bounds.at(0), domain_bounds.at(1)};
	double y_range[2] = {domain_bounds.at(2), domain_bounds.at(3)};
	double z_range[2] = {domain_bounds.at(4), domain_bounds.at(5)};
	Domain domain = Domain(x_range, y_range, z_range, contact_tol);
	double domain_vol = domain.volume();

	// Generates particles acording to the specification given in component-file
	std::vector<Superellipsoid*>* particles = generate_random_particles(components, domain_vol);

	// TODO: double check that the target volume fractions are achieved after the change ParticleDistribution -> Component
	// TODO: fill the domain using advancing front!!
	// TODO: Compute contact statistics and output it in some reasonable way
	// TODO: DONE!
	
	return EXIT_SUCCESS;
}
