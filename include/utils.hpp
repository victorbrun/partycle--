#ifndef UTILS_H
#define UTILS_H

#include "particle_generation.hpp"

#include <vector>
#include <string>

/**
 *
 */
std::vector<std::vector<std::string>> read_csv(const std::string& file_name);

/**
 *
 */
std::vector<std::string> string_split(const std::string& s, const std::string& delim);

/**
 * 
 */
std::vector<double> parse_domain(const std::string& domain_string);

/**
 *
 */
std::vector<ParticleDistribution> parse_components(const std::string& file_name);

/**
 *
 */

#endif
