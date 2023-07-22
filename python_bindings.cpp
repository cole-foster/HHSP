#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include "HierarchicalHSP.hpp"


class HHSP_Bindings {
public:
    // dataset info
    size_t dimension_{0};
    size_t dataset_size_{0};
    float* data_pointer_ = nullptr;

    // HHSP algo
    HierarchicalHSP* alg_hhsp{nullptr};

    // Constructors
    HHSP_Bindings(size_t dimension) {
        dimension_ = dimension;
        alg_hhsp = new HierarchicalHSP(dimension_);
        printf("--- hint: call `add_dataset()`\n");
    };
    ~HHSP_Bindings() {
        if (alg_hhsp) delete alg_hhsp;
    };

    inline size_t get_dimension() {return dimension_;}

    // Take Numpy Buffer as Input: allows for pointer to original dataset storage
    inline void add_dataset(py::buffer& input_data_buffer) {

        // process the incoming data
        py::array_t<float, py::array::c_style | py::array::forcecast> items(input_data_buffer);
        auto buffer = items.request();
        size_t rows, features;
        if (buffer.ndim == 2) {
            rows = buffer.shape[0];
            features = buffer.shape[1];
        } else std::runtime_error("Input array must be 2 dimensional");
        if (features != dimension_) std::runtime_error("Input Array Does Not Match Given Dimension");
        dataset_size_ = rows;
        printf("Given dataset of %u-dimensional vectors with N=%u\n",dimension_, dataset_size_);

        // add the data pointer to index, without duplicating memory
        data_pointer_ = static_cast<float*> (buffer.ptr);
        alg_hhsp->set_data_pointer(data_pointer_, dataset_size_);

        return;
    }

    inline void print_representation(elementID const member_index) {
        alg_hhsp->print_representation(member_index);
    }

    // Take Numpy Buffer as Input: allows for pointer to original dataset storage
    inline py::list HSP_Test_Members(elementID const member_index) {
        if (!data_pointer_) std::runtime_error("Dataset Not Initialized: Run `add_dataset()` before");
        if (member_index < 0 || member_index > dataset_size_) std::runtime_error("given index not in dataset");
        std::vector<elementID> neighbors{};
        alg_hhsp->HSP_Test(member_index, neighbors);
        py::list neighbors_py = py::cast(neighbors);
        return neighbors_py;
    }

    // Take Numpy Buffer as Input: allows for pointer to original dataset storage
    inline py::list HSP_Test(py::buffer& query_data_buffer) {
        // process the incoming data
        py::array_t<float, py::array::c_style | py::array::forcecast> items(query_data_buffer);
        auto buffer = items.request();
        size_t rows, features;
        if (buffer.ndim == 2) {
            rows = buffer.shape[0];
            features = buffer.shape[1];
        } else std::runtime_error("Input array must be 2 dimensional");
        if (features != dimension_) std::runtime_error("Input Array Does Not Match Given Dimension");
        int num_queries = rows;
        printf("Given testset of %u-dimensional vectors with %u queries\n",dimension_, num_queries);

        // Find HSP neighbors sequentially
        std::vector<std::vector<elementID>> hsp_neighbors{};
        hsp_neighbors.resize(num_queries);
        for (int i = 0; i < num_queries; i++) {
            alg_hhsp->HSP_Test((float*)items.data(i), hsp_neighbors[i]);
        }

        // ParallelFor(0, rows, num_threads, [&](size_t row, size_t threadId) {
        // std::priority_queue<std::pair<dist_t, hnswlib::labeltype>> result =
        //     appr_alg->searchKnn((void*)items.data(row), k, p_idFilter);
        // if (result.size() != k)
        //     throw std::runtime_error(
        //         "Cannot return the results in a contigious 2D array. Probably ef or M is too small");
        // for (int i = k - 1; i >= 0; i--) {
        //     auto& result_tuple = result.top();
        //     data_numpy_d[row * k + i] = result_tuple.first;
        //     data_numpy_l[row * k + i] = result_tuple.second;
        //     result.pop();
        // }
        py::list out = py::cast(hsp_neighbors);
        return out;
    }
};


PYBIND11_MODULE(HHSP, handle) {
    handle.doc() = "Python Module for Hierarchical HSP";

    py::class_<HHSP_Bindings> (handle, "HHSP")
        .def(py::init<size_t>())
        .def_property_readonly("dimension", &HHSP_Bindings::get_dimension)
        .def("add_dataset", &HHSP_Bindings::add_dataset, py::arg("data buffer"))
        .def("print_representation", &HHSP_Bindings::print_representation, py::arg("element"))
        .def("HSP_Test_Members", &HHSP_Bindings::HSP_Test_Members, py::arg("member index"))
        .def("HSP_Test", &HHSP_Bindings::HSP_Test, py::arg("query data"));
        // .def("multiply", &HHSP_Bindings::multiply);
}