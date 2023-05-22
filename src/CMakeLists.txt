# Sparse-Matrix: for distance computations
add_library(sparse_matrix sparse-matrix.cpp sparse-matrix.hpp)
target_include_directories(sparse_matrix PUBLIC ${PROJECT_SOURCE_DIR}/include/sparse-map/include)

# Pivot-Layer: for constructing the index
add_library(pivot_layer pivot-layer.hpp pivot-index.hpp pivot-index.cpp)
target_link_libraries(pivot_layer PUBLIC sparse_matrix)

# GHSP: for searching
add_library(src GHSP.hpp GHSP.cpp NNS.hpp NNS.cpp)
target_link_libraries(src PUBLIC pivot_layer)
