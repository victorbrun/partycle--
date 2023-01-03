#######################################################
### Partycle--                                       ###
#######################################################
# Project information
cmake_minimum_required(VERSION 3.15)
project(Partycle-- VERSION 0.0.1)

# Includes
include(FetchContent)

# Fetches Matplot++
FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master) # or whatever tag you want

FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()


# Compiler setting 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

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

target_link_libraries(partycle-- PUBLIC matplot)
