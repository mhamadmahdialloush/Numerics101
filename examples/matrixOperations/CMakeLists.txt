cmake_minimum_required(VERSION 3.10)

project(Numerics101 VERSION 1.0)

# Set C++ standard (optional but recommended)
set(CMAKE_CXX_STANDARD 20)

# Find MPI package
find_package(MPI REQUIRED)

# Set the executable names
set(example1 example1.exe)
set(example2 example2.exe)
set(example3 example3.exe)

# Find all source files in the src directory to a variable
file(GLOB_RECURSE SOURCES "../../src/*.cpp")

# Find all source files in the src directory to a variable
file(GLOB_RECURSE HEADER_FILES "../../src/*.h" "../../src/*.hpp")

# Add the executable target for example1
add_executable(${example1} example1.cpp ${SOURCES})

# Add the executable target example2
add_executable(${example2} example2.cpp ${SOURCES})

# Add the executable target example3
add_executable(${example3} example3.cpp ${SOURCES})

# Include the directories containing the header files
foreach(header_file ${HEADER_FILES})
    get_filename_component(header_dir ${header_file} DIRECTORY)
    target_include_directories(${example1} PRIVATE ${header_dir})
    target_include_directories(${example2} PRIVATE ${header_dir})
    target_include_directories(${example3} PRIVATE ${header_dir})
endforeach()


# Add MPI include directories to your project
include_directories(${MPI_INCLUDE_PATH})

# Link ser executable with MPI libraries
target_link_libraries(${example1} PRIVATE ${MPI_CXX_LIBRARIES})

# Link par executable with MPI libraries
target_link_libraries(${example2} PRIVATE  ${MPI_CXX_LIBRARIES})

# Link par executable with MPI libraries
target_link_libraries(${example3} PRIVATE  ${MPI_CXX_LIBRARIES})

