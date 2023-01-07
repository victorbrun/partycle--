#include "program_options.hpp"
#include <stdexcept>
#include <vector>
#include <iostream>

namespace {
	static std::vector<std::string> _input_args;
}

void program_options::parse(int argc, char* argv[]) {
	if (argc > 64) {
		throw std::runtime_error("to many input parameters.");
	}

	const std::vector<std::string> args(argv + 1, argv + argc);

	for (const std::string& arg : args) {
		_input_args.push_back(arg);		
	}
}

std::string program_options::get_option(const std::string& option_name) {
	for (auto it = _input_args.begin(), end = _input_args.end(); it != end; it++) {
		if (*it == option_name) {
			if (it + 1 != end) {
				return *(it + 1);
			}
		}
	}

	return "";
}
