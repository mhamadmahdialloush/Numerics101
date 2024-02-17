cmake_minimum_required(VERSION 3.10)

project(Numerics101 VERSION 1.0)

# Set C++ standard (optional but recommended)
set(CMAKE_CXX_STANDARD 20)

# Find MPI package
find_package(MPI REQUIRED)

# Set the executable names
set(ser_exec ser_petsc_example.exe)
set(par_exec par_petsc_example.exe)

# PkgConfig
find_package(PkgConfig)

# PETSc
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
endif()

if(PETSC_FOUND)
    add_definitions(-DUSE_PETSC)  # Define a preprocessor macro
    include_directories(${PETSC_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Error: PETSc is not found on your system.")
endif()

# Find all source files in the src directory to a variable
file(GLOB_RECURSE SOURCES "../../src/*.cpp")

# Find all source files in the src directory to a variable
file(GLOB_RECURSE HEADER_FILES "../../src/*.h" "../../src/*.hpp")

# Add the executable target for the serial exec
add_executable(${ser_exec} petscSolveSerBlockCSR.cpp ${SOURCES})

# Add the executable target fort the par exec
add_executable(${par_exec} petscSolveParBlockCSR.cpp ${SOURCES})


# Include the directories containing the header files
foreach(header_file ${HEADER_FILES})
    get_filename_component(header_dir ${header_file} DIRECTORY)
    target_include_directories(${ser_exec} PRIVATE ${header_dir})
    target_include_directories(${par_exec} PRIVATE ${header_dir})
endforeach()


# Add MPI include directories to your project
include_directories(${MPI_INCLUDE_PATH})

# Link ser executable with MPI libraries
target_link_libraries(${ser_exec} PRIVATE ${MPI_CXX_LIBRARIES})

# Link par executable with MPI libraries
target_link_libraries(${par_exec} PRIVATE  ${MPI_CXX_LIBRARIES})

# Link to PETSc libraries
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(${ser_exec} PRIVATE ${PETSC_INCLUDE_DIRS}/../lib/libpetsc.3.20.1.dylib)
    target_link_libraries(${par_exec} PRIVATE ${PETSC_INCLUDE_DIRS}/../lib/libpetsc.3.20.1.dylib)
ELSE()
    target_link_libraries(${ser_exec} PRIVATE ${PETSC_LIBRARIES})
    target_link_libraries(${par_exec} PRIVATE ${PETSC_LIBRARIES})
ENDIF()

