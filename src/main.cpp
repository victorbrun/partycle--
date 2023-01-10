#include "particle_generation.hpp"
#include "program_options.hpp"
#include "superellipsoid.hpp"
#include <cctype>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <filesystem>
#include <string>

std::vector<std::string> string_split(const std::string& s, const std::string& delim) {
	// Allocates output vector and finds first 
	// occurance of delimitier too prepair for while loop
	std::vector<std::string> vec;
	std::string sub_string = s;
	int delim_idx = s.find(delim);

	while (delim_idx > -1) {
		// Adds the string up until the delimiter to output vector
		std::string split_string = sub_string.substr(0, delim_idx);
		vec.push_back(split_string);

		// Updates sub_string to remove the sub string added to vec,
		// and find the index for the occurance of the next delimiter.
		sub_string = sub_string.substr(delim_idx + 1, sub_string.size());	
		delim_idx = sub_string.find(delim);
	}

	// Above loop will leave the last part of the string in sub_string.
	// If this is not an empty string we want to add it to the output vector.
	if (sub_string != "") vec.push_back(sub_string);

	return vec;
}

std::vector<double> parse_domain(const std::string& domain_string) {
	// Checks that the domain definition has the right format
	std::regex domain_string_format("^\\[[0-9]+\\,[0-9]+\\]x\\[[0-9]+\\,[0-9]+\\]x\\[[0-9]+\\,[0-9]+\\]$");
	if (!std::regex_match(domain_string, domain_string_format)) {
		throw std::runtime_error("particle--: invalid domain definition format. It should have the format: [ax,bx]x[ay,by]x[az,bz], where ai is lower bound for i-axis and bi is upper bound for the same axis.");
	}
	
	// Allocates vector in which to store bounds
	std::vector<double> domain_bounds(6);
	
	// The result after splitting will be {"[ax, bx]", "[ay,by]", "[az,bz]"}.
	// So by removing first and last element in exh string an splitting by ","
	// we extract the bounds.
	std::vector<std::string> split1_vec = string_split(domain_string, "x");
	for (size_t ix = 0; ix < split1_vec.size(); ix++) {
		std::string temp_string = split1_vec.at(ix).substr(1, split1_vec.at(ix).size()-1);
		std::vector<std::string> split2_vec = string_split(temp_string, ",");

		// Converting string to double and assigning it to vector entries
		domain_bounds.at(2*ix + 0) = std::stod(split2_vec.at(0));
		domain_bounds.at(2*ix + 1) = std::stod(split2_vec.at(1));
	}

	return domain_bounds;
}

std::vector<ParticleDistribution> parse_components_file(const std::string& file_name) {
	// Checking that the component-file is a csv file 
	std::regex file_format1(".*\\.csv$");
	std::regex file_format2(".*\\.CSV$");
	if ( !std::regex_match(file_name, file_format1) && !std::regex_match(file_name, file_format2) ) {
		throw std::runtime_error("particle--: specified component-file must be csv file without spaces in name, i.e. have .csv or .CSV file extension.");
	}

	// Checks if component-file exists
	if (!std::filesystem::exists(file_name)) {
		std::runtime_error("particle--: file does not exist.");
	}

	std::vector<std::vector<std::string>> content = read_csv(file_name);

	std::vector<ParticleDistribution> pds;
	return pds;

}

std::vector<std::vector<std::string>> read_csv(const std::string& file_name) {
	std::fstream file(file_name, std::ios::in);
	if (!file.is_open()) {
		throw std::runtime_error("particle--: could not open file.");
	}
	
	std::vector<std::vector<std::string>> rows;
	std::string line;
	size_t ix = 0;
	while (std::getline(file, line)) {
		// Extracting each field for the given line
		std::vector<std::string> entries = string_split(line, ";");
		for (size_t jx = 0; jx < entries.size(); jx++) {
			// Removes whitespace of each entry
			std::regex white_space("\\s+");
			entries.at(ix) = std::regex_replace(entries.at(ix), white_space, "");
			std::cout << entries.at(ix) << std::endl;
			rows.push_back(entries);
		}
	}

	return rows;
}

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
	std::vector<ParticleDistribution> pds;
	try {
		pds = parse_components_file(components_file_name);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << "\n";
		return EXIT_FAILURE;
	}

	// TODO: read the file and parse the contents
	// TODO: initiate domain and generate particles
	// TODO: fill the domain using advancing front!!
	// TODO: Compute contact statistics and output it in some reasonable way
	// TODO: DONE!
	
	return EXIT_SUCCESS;
}
