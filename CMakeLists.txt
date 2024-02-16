cmake_minimum_required(VERSION 3.10)

project(Numerics101 VERSION 1.0)

# Set C++ standard (optional but recommended)
set(CMAKE_CXX_STANDARD 20)

# Find MPI package
find_package(MPI REQUIRED)

# Set the executable name
set(exec numerics_exec.exe)

# PkgConfig
find_package(PkgConfig)

# PETSc
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
endif()

if(PETSC_FOUND)
    message("PETSc functionalities are enabled.")
    add_definitions(-DUSE_PETSC)  # Define a preprocessor macro
    include_directories(${PETSC_INCLUDE_DIRS})
else()
    # Include PETSc include directories
    if(DEFINED PETSC_DIR)
        message("PETSC_DIR (${PETSC_DIR}) is provided. PETSc functionalities are enabled.")
        add_definitions(-DUSE_PETSC)  # Define a preprocessor macro
        include_directories(${PETSC_DIR}/include)

        if(DEFINED PETSC_ARCH)
            include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)
        endif()
    else()
        message("PETSc functionalities are disabled.")
    endif()
endif()

# Add all source files in the source directory to a variable
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Add the executable target
add_executable(${exec} main.cpp ${SOURCES})

# Add MPI include directories to your project
include_directories(${MPI_INCLUDE_PATH})

# Specify the include directories for the target
target_include_directories(${exec} PRIVATE include)

# Link your executable with MPI libraries
target_link_libraries(${exec} ${MPI_CXX_LIBRARIES})

# Link to PETSc libraries
if(PETSC_FOUND)
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
        target_link_libraries(${exec} PRIVATE ${PETSC_INCLUDE_DIRS}/../lib/libpetsc.3.20.1.dylib)
    ELSE()
        target_link_libraries(${exec} ${PETSC_LIBRARIES})
    ENDIF()
else()
    if(DEFINED PETSC_DIR)
        if(DEFINED PETSC_ARCH)
            target_link_libraries(${exec} PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)
        else()
            target_link_libraries(${exec} PRIVATE ${PETSC_DIR}/lib/libpetsc.so)
        endif()
    endif()
endif()