cmake_minimum_required(VERSION 3.10)

project(Numerics101 VERSION 1.0)

# Set the executable name
set(exec numerics_exec.exe)

# Add all source files in the source directory to a variable
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Add the executable target
add_executable(${exec} main.cpp ${SOURCES})

# Specify the include directories for the target
target_include_directories(${exec} PRIVATE include)

# Set the C++ standard for the target
target_compile_features(${exec} PRIVATE cxx_std_17)
