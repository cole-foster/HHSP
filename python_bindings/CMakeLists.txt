# add the pybind files
add_subdirectory(pybind11)
pybind11_add_module(HHSP python_bindings.cpp)
target_link_libraries(HHSP PUBLIC SRC)
target_include_directories(HHSP PUBLIC ${PROJECT_SOURCE_DIR}/src)
