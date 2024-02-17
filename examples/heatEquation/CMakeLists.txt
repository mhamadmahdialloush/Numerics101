cmake_minimum_required(VERSION 3.10)

project(Numerics101 VERSION 1.0)

# Set C++ standard (optional but recommended)
set(CMAKE_CXX_STANDARD 20)

# Find MPI package
find_package(MPI REQUIRED)

# Set the executable names
set(exec solveHeatEquation.exe)

# Find all source files in the src directory to a variable
file(GLOB_RECURSE SOURCES "../../src/*.cpp")

# Find all source files in the src directory to a variable
file(GLOB_RECURSE HEADER_FILES "../../src/*.h" "../../src/*.hpp")

# Add the executable target for exec
add_executable(${exec} solveHeatEquation.cpp ${SOURCES})

# Include the directories containing the header files
foreach(header_file ${HEADER_FILES})
    get_filename_component(header_dir ${header_file} DIRECTORY)
    target_include_directories(${exec} PRIVATE ${header_dir})
endforeach()


# Add MPI include directories to your project
include_directories(${MPI_INCLUDE_PATH})

# Link ser executable with MPI libraries
target_link_libraries(${exec} PRIVATE ${MPI_CXX_LIBRARIES})

