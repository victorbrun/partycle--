#include "program_options.hpp"
#include <cstdlib>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <filesystem>

int main(int argc, char* argv[]) {
	// Parses input arguments and stores them in file-scope varaibel	
	program_options::parse(argc, argv);

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

	// Checks that the domain definition has the right format
	std::regex domain_string_format("^\\[[0-9]+\\,[0-9]+\\]x\\[[0-9]+\\,[0-9]+\\]x\\[[0-9]+\\,[0-9]+\\]$");
	if (!std::regex_match(domain_string, domain_string_format)) {
		std::cerr << 
			"particle--: invalid domain definition format. It should have the format: [ax,bx]x[ay,by]x[az,bz], where ai is lower bound for i-axis and bi is upper bound for the same axis." <<
			"\n";
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

	// Checking that the component-file is a csv file 
	
	std::regex file_format1("\\w*\\.csv$");
	std::regex file_format2("\\w*\\.CSV$");
	if ( !std::regex_match(components_file_name, file_format1) && !std::regex_match(components_file_name, file_format2) ) {
		std::cerr << "particle--: specified component-file must be csv file without spaces in name, i.e. have .csv or .CSV file extension." << "\n";
		return EXIT_FAILURE;
	}

	// Checks if component-file exists
	if (!std::filesystem::exists(components_file_name)) {
		std::cerr << "particle--: the file: " << components_file_name << " does not exist." << "\n";
		return EXIT_FAILURE;
	}

	// TODO: parse the domain bounds 
	// TODO: read the file and parse the contents
	// TODO: initiate domain and generate particles
	// TODO: fill the domain using advancing front!!
	// TODO: Compute contact statistics and output it in some reasonable way
	// TODO: DONE!
	
	return EXIT_SUCCESS;
}
