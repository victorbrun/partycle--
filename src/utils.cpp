#include "utils.hpp"
#include "particle_generation.hpp"
#include "superellipsoid.hpp"
#include <fstream>
#include <regex>
#include <iostream>
#include <filesystem>

std::vector<std::vector<std::string>> read_csv(const std::string& file_name) {
	std::fstream file(file_name, std::ios::in);
	if (!file.is_open()) {
		throw std::runtime_error("particle--: could not open file.");
	}
	
	std::vector<std::vector<std::string>> rows;
	std::string line;
	while (std::getline(file, line)) {
		// Extracting each field for the given line
		std::vector<std::string> entries = string_split(line, ";");
		rows.push_back(entries);
	}

	return rows;
}

std::vector<std::string> string_split(const std::string& s, const std::string& delim) {
	// This function would probably be more efficient by streaming the 
	// the string, but I do not know how that works and I cannot bother 
	// to learn that right now :)

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


std::vector<ParticleDistribution> parse_components(const std::string& file_name) {
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

	// Reads csv into nested vector of strings and checks that the first row is 
	// a correctly formated header
	std::vector<std::vector<std::string>> content = read_csv(file_name);
	if (content.size() < 2) {
		throw std::runtime_error("particle--: component-file must have at least two rows where the first one is header.");
	} else if (content.at(0).size() != 7) {
		throw std::runtime_error("particle--: component-file must have the header: class;volume_distribution;reference_particle_a;reference_particle_b;reference_particle_c;reference_particle_n1;reference_particle_n2");
	}

	// Finally checks that the header has the right exact format with spelling and column name placement
	std::vector<std::string> header = content.at(0);
	bool header_check = (header.at(0) == "class") && (header.at(1) == "volume_distribution") && (header.at(2) == "reference_particle_a") && 
						(header.at(3) == "reference_particle_b") && (header.at(4) == "reference_particle_c") && (header.at(5) == "reference_particle_n1") && (header.at(6) == "reference_particle_n2");
	if (!header_check) {
		throw std::runtime_error("particle--: component-file must have the header: class;volume_distribution;reference_particle_a;reference_particle_b;reference_particle_c;reference_particle_n1;reference_particle_n2");
	}	
	
	// Extracting the information from contents and turning it into ParticleDistribution
	std::vector<ParticleDistribution> pds; // For some reason i could not pre-allocate the size of pds so we push onto it
	for (size_t ix = 1; ix < content.size(); ix++) {
		std::vector<std::string> row = content.at(ix);
		int cls = std::stoi(row.at(0));
		Distribution vol_dist = parse_distribution(row.at(1));
		double scale_params[3] = {std::stod(row.at(2)), std::stod(row.at(3)), std::stod(row.at(4))};
		double shape_params[2] = {std::stod(row.at(5)), std::stod(row.at(6))};
			
		Superellipsoid p = Superellipsoid(cls, scale_params, shape_params);
		ParticleDistribution pd = {cls, p, vol_dist};
		pds.push_back(pd);
	}

	return pds;
}

