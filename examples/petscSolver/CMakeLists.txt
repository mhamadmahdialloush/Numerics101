cmake_minimum_required(VERSION 3.10)

project(Numerics101 VERSION 1.0)

# Set C++ standard (optional but recommended)
set(CMAKE_CXX_STANDARD 20)

# Find MPI package
find_package(MPI REQUIRED)

# Set the executable names
set(ser_exec ser_petsc_example.exe)
set(ser_block_exec ser_petsc_block_example.exe)
set(par_block_exec par_petsc_block_example.exe)

# PkgConfig
find_package(PkgConfig)

# Find all source files in the src directory to a variable
file(GLOB_RECURSE SOURCES "../../src/*.cpp")

# Find all source files in the src directory to a variable
file(GLOB_RECURSE HEADER_FILES "../../src/*.h" "../../src/*.hpp")

# Add the executable target for the serial exec
add_executable(${ser_exec} petscSolveSerCSR.cpp ${SOURCES})

# Add the executable target for the serial exec
add_executable(${ser_block_exec} petscSolveSerBlockCSR.cpp ${SOURCES})

# Add the executable target fort the par exec
add_executable(${par_block_exec} petscSolveParBlockCSR.cpp ${SOURCES})


# Include the directories containing the header files
foreach(header_file ${HEADER_FILES})
    get_filename_component(header_dir ${header_file} DIRECTORY)
    include_directories(${header_dir})
endforeach()

# PETSc
if(DEFINED PETSC_DIR)
    message("PETSC_DIR (${PETSC_DIR}) is provided.")
    add_definitions(-DUSE_PETSC)  # Define a preprocessor macro
    target_include_directories(${ser_exec} PRIVATE ${PETSC_DIR}/include)
    target_include_directories(${ser_block_exec} PRIVATE ${PETSC_DIR}/include)
    target_include_directories(${par_block_exec} PRIVATE ${PETSC_DIR}/include)

    if(DEFINED PETSC_ARCH)
        target_include_directories(${ser_exec} PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/include)
        target_include_directories(${ser_block_exec} PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/include)
        target_include_directories(${par_block_exec} PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/include)
    endif()
else()
    if (PKG_CONFIG_FOUND)
        pkg_check_modules(PETSC PETSc)

        if(PETSC_FOUND)
        target_include_directories(${ser_exec} PRIVATE ${PETSC_INCLUDE_DIRS})
	    target_include_directories(${ser_block_exec} PRIVATE ${PETSC_INCLUDE_DIRS})
	    target_include_directories(${par_block_exec} PRIVATE ${PETSC_INCLUDE_DIRS})
            add_definitions(-DUSE_PETSC)  # Define a preprocessor macro
        else()
	    message(FATAL_ERROR "Error: PETSc is not found on your system.")
        endif()
    endif()
endif()

# Add MPI include directories to your project
include_directories(${MPI_INCLUDE_PATH})

# Link ser executable with MPI libraries
target_link_libraries(${ser_exec} PRIVATE ${MPI_CXX_LIBRARIES})

# Link ser executable with MPI libraries
target_link_libraries(${ser_block_exec} PRIVATE ${MPI_CXX_LIBRARIES})

# Link par executable with MPI libraries
target_link_libraries(${par_block_exec} PRIVATE ${MPI_CXX_LIBRARIES})

# Link to PETSc libraries
if(DEFINED PETSC_DIR)
    if(DEFINED PETSC_ARCH)
        target_link_libraries(${ser_exec} PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)
        target_link_libraries(${ser_block_exec} PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)
        target_link_libraries(${par_block_exec} PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)
    else()
        target_link_libraries(${ser_exec} PRIVATE ${PETSC_DIR}/lib/libpetsc.so)
        target_link_libraries(${ser_block_exec} PRIVATE ${PETSC_DIR}/lib/libpetsc.so)
        target_link_libraries(${par_block_exec} PRIVATE ${PETSC_DIR}/lib/libpetsc.so)
    endif()
else()
    if(PETSC_FOUND)
        IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
            target_link_libraries(${ser_exec} PRIVATE ${PETSC_INCLUDE_DIRS}/../lib/libpetsc.3.20.1.dylib)
            target_link_libraries(${ser_block_exec} PRIVATE ${PETSC_INCLUDE_DIRS}/../lib/libpetsc.3.20.1.dylib)
            target_link_libraries(${par_block_exec} PRIVATE ${PETSC_INCLUDE_DIRS}/../lib/libpetsc.3.20.1.dylib)
        ELSE()
            target_link_libraries(${ser_exec} PRIVATE ${PETSC_LIBRARIES})
            target_link_libraries(${ser_block_exec} PRIVATE ${PETSC_LIBRARIES})
            target_link_libraries(${par_block_exec} PRIVATE ${PETSC_LIBRARIES})
        ENDIF()
    endif()
endif()
