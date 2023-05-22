cmake_minimum_required(VERSION 3.15)
project(GHSP
        VERSION 1.0
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# User defined options
option(CMAKE_BUILD_TYPE "type of build [Release/Debug/RelWithDebInfo]" Release)

# -- Crucial flags for performance --
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native" CACHE STRING "" FORCE)

# -- This one is the most important one - order of magnitude faster --
# -- See this: https://gcc.gnu.org/ml/gcc/2001-07/msg02150.html --
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" CACHE STRING "" FORCE)

# openmp
find_package(OpenMP)
MESSAGE("OpenMP: ${OpenMP_pthread_LIBRARY}")

# GHSP optimization
option(GHSP_OPT "enable use/build of optimization library" OFF)
if (GHSP_OPT)
    SET(NLOPT_HOME "/users/cfoste18/data/cfoste18/hRNG/hRNG_Optimization/include/nlopt/install")
    SET(NLOPT_DIRECTORIES "${NLOPT_HOME}/include/")
    find_library(NLOPT nlopt PATHS "${NLOPT_HOME}/lib/")
endif()
MESSAGE("NLOPT: ${NLOPT}")

# set CMAKE_INSTALL_PREFIX as default location of `make install` command
# creates a bin/ folder within the location
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# add the source files
add_subdirectory(src)

# add the scripts
add_subdirectory(scripts)



