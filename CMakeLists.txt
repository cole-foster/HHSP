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

# hRNG optimization
set(HRNG_OPT TRUE)
if (HRNG_OPT)
    #SET(NLOPT_HOME "/home/colef8/nlopt")
    SET(NLOPT_HOME "/users/cfoste18/data/cfoste18/hRNG/hRNG_Optimization/include/nlopt/install")
    SET(NLOPT_DIRECTORIES "${NLOPT_HOME}/include/")
    find_library(NLOPT nlopt PATHS "${NLOPT_HOME}/lib/")
endif()
MESSAGE("NLOPT: ${NLOPT}")


# set CMAKE_INSTALL_PREFIX as default location of `make install` command
# creates a bin/ folder within the location
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# all files
set(HEADERS CLI11.hpp datasets.hpp sparse-matrix.hpp pivot.hpp pivot-index.hpp NNS.hpp GHSP.hpp)
set(SOURCES sparse-matrix.cpp pivot-index.cpp NNS.cpp GHSP.cpp)

# add scripts
add_executable(main main.cpp ${HEADERS} ${SOURCES})
install(TARGETS main)

add_executable(opt main-optimization.cpp ${HEADERS} ${SOURCES})
install(TARGETS opt)
target_include_directories(opt PRIVATE ${NLOPT_DIRECTORIES})
target_link_libraries(opt PRIVATE ${NLOPT})
