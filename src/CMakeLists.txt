# Sparse-Matrix: for distance computations
# add_library(sparse_matrix sparse-matrix.cpp sparse-matrix.hpp)
# target_include_directories(sparse_matrix PUBLIC ${PROJECT_SOURCE_DIR}/include/sparse-map/include)

# Pivot-Layer: for constructing the index

# # GHSP: for searching
# add_library(src GHSP.hpp GHSP.cpp NNS.hpp NNS.cpp)
# target_link_libraries(src PUBLIC pivot_layer)

# HHSP: overall library
set(SRC_SOURCES HierarchicalHSP.hpp HierarchicalHSP.cpp metrics.hpp pivot-layer.hpp)
add_library(src SHARED ${SRC_SOURCES})
target_include_directories(src PUBLIC ${PROJECT_SOURCE_DIR}/include/sparse-map/include)
set_target_properties(src PROPERTIES POSITION_INDEPENDENT_CODE ON) # not sure why... needed for shared libs?
            # https://stackoverflow.com/questions/19364969/compilation-fails-with-relocation-r-x86-64-32-against-rodata-str1-8-can-not