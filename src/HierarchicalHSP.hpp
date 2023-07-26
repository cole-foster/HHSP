#ifndef HierarchicalHSP_hpp
#define HierarchicalHSP_hpp

#include <algorithm>
#include <cstdio>
#include <memory>
#include <vector>
#include <tsl/sparse_map.h>

#include "metrics.hpp"
#include "pivot-layer.hpp"

/*
    Heavily influenced by https://github.com/nmslib/hnswlib
*/
typedef unsigned int elementID;  // elementID is 4-bytes always

/**
 * @brief Hierarchical HSP
 *
 */
class HierarchicalHSP {
   private:
    int num_threads_ = 1;

    // dataset storage
    size_t dimension_ = 0;
    size_t dataset_size_ = 0;
    float* data_pointer_ = nullptr;

    //> get the pointer to the beginning of index's representation in data_pointer_
    inline float* internalDataPointer_(elementID const index) const { return data_pointer_ + index * dimension_; }

    // distance metric: from metrics.h, needs void*
    DISTFUNC<float> distance_function_;
    void* dist_func_dimension_ = nullptr;

    // index
    std::shared_ptr<std::vector<PivotLayer>> index_ = std::make_shared<std::vector<PivotLayer>>();

   public:
    HierarchicalHSP(){};
    HierarchicalHSP(size_t const dimension) : dimension_(dimension) {
        distance_function_ = metrics::EuclideanDistance;
        dist_func_dimension_ = &dimension_;
    };
    HierarchicalHSP(size_t const dimension, size_t const dataset_size, float*& data_pointer)
        : dimension_(dimension), dataset_size_(dataset_size) {
        distance_function_ = metrics::EuclideanDistance;
        dist_func_dimension_ = &dimension_;
        data_pointer_ = data_pointer;  // char: 1 byte per element, float: 4 bytes per el
    };
    ~HierarchicalHSP() {};

    inline void set_data_pointer(float*& data_pointer, size_t const dataset_size) {
        data_pointer_ = data_pointer;
        dataset_size_ = dataset_size;
    }

    //> compute the distance between two points
    inline float compute_distance(elementID const index1, elementID const index2) const {
        return distance_function_(internalDataPointer_(index1), internalDataPointer_(index2), dist_func_dimension_);
    }
    inline float compute_distance(float* index1_ptr, float* index2_ptr) const {
        return distance_function_(index1_ptr, index2_ptr, dist_func_dimension_);
    }

    //> perform the actual HSP Test: can be done in parallel
    void HSP_Test(float* query_pointer, std::vector<elementID>& neighbors) const;
    inline void HSP_Test(elementID const index, std::vector<elementID>& neighbors) const {
        HSP_Test(internalDataPointer_(index), neighbors);
    }

    //> Creating the Index for HHSP
    void create_index(std::vector<float> radius_vector);

    //> Performing HHSP Search
    void Hierarchical_HSP_Test(float* query_pointer, std::vector<elementID>& neighbors) const;
    inline void Hierarchical_HSP_Test(elementID const index, std::vector<elementID>& neighbors) const {
        Hierarchical_HSP_Test(internalDataPointer_(index), neighbors);
    }
    float get_distance(elementID const index1, float* index1_ptr, elementID const index2,
                       tsl::sparse_map<int, std::vector<float>>& distance_store_array) const;
    bool validate_point(elementID const point_index, elementID const query_index, float* query_ptr,
                        std::vector<elementID> const& neighbors,
                        tsl::sparse_map<int, std::vector<float>>& distance_store_array) const;
    bool validate_pivot(elementID const pivot_index, float const radius, elementID const query_index, float* query_ptr,
                        std::vector<elementID> const& neighbors,
                        tsl::sparse_map<int, std::vector<float>>& distance_store_array, int& pivot_outcome) const;

    //          OUTPUT AND DEBUGGING

    //> print the entire representation for the index
    inline void print_representation(elementID const index) {
        float* index_pointer = internalDataPointer_(index);
        printf("%u: [", (unsigned int)index);
        for (int i = 0; i < (int)dimension_; i++) {
            printf("%.6f,", *index_pointer);
            index_pointer++;
        }
        printf("]\n");
        return;
    }

    inline void print_vector(std::vector<elementID>& vec) const {
        printf("[");
        for (int i = 0; i < (int)vec.size(); i++) {
            printf("%u, ", vec[i]);
        }
        printf("]\n");
        return;
    }

    inline void print_hierarchy_stats() const {
        int num_layers = (int)index_->size() + 1;
        printf("%d-Layer Hierarchy on N=%u in %d-D\n", num_layers, (unsigned int)dataset_size_, (int)dimension_);
        int layer_index = 0;
        for (; layer_index < num_layers - 1; layer_index++) {
            unsigned int num_pivots = (*index_)[layer_index].get_num_pivots();
            float radius = (*index_)[layer_index].get_radius();
            printf("  %d: %.4f --> %u\n", layer_index, radius, num_pivots);
        }
        printf("  %d: %.4f --> %u\n", layer_index, 0.00, (unsigned int)dataset_size_);
    }

    // validating the hierarchy: ensures minimal coverage
    bool validate_index();
    bool recursive_depth_first_check(std::vector<unsigned int> parents, int const layer_index,
                                     std::vector<unsigned int>& points_list);
};

#endif  // HierarchicalHSP_hpp