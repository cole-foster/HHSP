# HHSP: overall library
set(SRC_SOURCES HierarchicalHSP.hpp HierarchicalHSP.cpp metrics.hpp pivot-layer.hpp)
add_library(SRC SHARED ${SRC_SOURCES})
target_include_directories(SRC PUBLIC ${PROJECT_SOURCE_DIR}/include/sparse-map/include)
set_target_properties(SRC PROPERTIES POSITION_INDEPENDENT_CODE ON) # not sure why... needed for shared libs?
            # https://stackoverflow.com/questions/19364969/compilation-fails-with-relocation-r-x86-64-32-against-rodata-str1-8-can-not