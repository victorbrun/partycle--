#ifndef PROGRAM_OPTIONS_H
#define PROGRAM_OPTIONS_H

#include <vector>
#include <string>

namespace program_options {
	void parse(int argc, char* argv[]);

	std::string get_option(const std::string& option_name);
}

#endif
