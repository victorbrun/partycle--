#######################################################
### Partycle--                                       ###
#######################################################
# Project information
cmake_minimum_required(VERSION 3.15)
project(Partycle-- VERSION 0.0.1)

# Compiler setting 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")

# Binary location
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()

# Setting up executable
set(SOURCES 
	src/main.cpp 

	src/program_options.cpp
	src/program_options.hpp

	src/utils.cpp
	src/utils.hpp

	src/superellipsoid.cpp  
	src/superellipsoid.hpp 
		
	src/coordinate_indexer.cpp 
	src/coordinate_indexer.hpp 
			
	src/domain.cpp 
	src/domain.hpp 
			
	src/particle_generation.cpp 
	src/particle_generation.hpp
)
add_executable(partycle-- ${SOURCES})
