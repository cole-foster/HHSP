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

# set CMAKE_INSTALL_PREFIX as default location of `make install` command
# creates a bin/ folder within the location
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# all files
set(HEADERS CLI11.hpp datasets.hpp sparse-matrix.hpp pivot-layer.hpp pivot-index.hpp NNS.hpp GHSP.hpp)
set(SOURCES sparse-matrix.cpp pivot-index.cpp NNS.cpp GHSP.cpp)
add_library(src ${HEADERS} ${SOURCES})
target_include_directories(src PUBLIC ${PROJECT_SOURCE_DIR}/sparse-map/include)

# add scripts
add_executable(main main.cpp)
target_link_libraries(main src)
install(TARGETS main)

# add_executable(opt main-optimization.cpp ${HEADERS} ${SOURCES})
# target_include_directories(opt PUBLIC ${PROJECT_SOURCE_DIR}/sparse-map/include)
# install(TARGETS opt)