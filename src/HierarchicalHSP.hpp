#ifndef HierarchicalHSP_hpp
#define HierarchicalHSP_hpp

#include <cstdio>
#include <vector>
#include "metrics.h"

/*
    Heavily influenced by https://github.com/nmslib/hnswlib
*/
typedef unsigned int elementID; // elementID is 4-bytes always

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
        // char* data_pointer_ = nullptr;              // char for indexing by bytes
        float* data_pointer_ = nullptr;
        // size_t size_data_per_element_ = 0;          // 

        //> get the pointer to the beginning of index's representation in data_pointer_
        inline float* internalDataPointer_(elementID const index) const {
            return data_pointer_ + index * dimension_;
        }
        
        // distance metric: from metrics.h, needs void*
        DISTFUNC<float> distance_function_;
        void *dist_func_dimension_ = nullptr;

    public:
        HierarchicalHSP() {};
        HierarchicalHSP(size_t const dimension): dimension_(dimension) {
            distance_function_ = metrics::EuclideanDistance;
            dist_func_dimension_ = &dimension_;
        };
        HierarchicalHSP(size_t const dimension, size_t const dataset_size, float*& data_pointer): dimension_(dimension),
                dataset_size_(dataset_size) {
            distance_function_ = metrics::EuclideanDistance;
            dist_func_dimension_ = &dimension_;
            data_pointer_ = data_pointer; // char: 1 byte per element, float: 4 bytes per el
        };
        ~HierarchicalHSP(){};

        //> compute the distance between two points
        inline float compute_distance(elementID const index1, elementID const index2) const {
            return distance_function_(internalDataPointer_(index1), internalDataPointer_(index2), dist_func_dimension_);
        }
        inline float compute_distance(float* index1_ptr, float* index2_ptr) const {
            return distance_function_(index1_ptr, index2_ptr, dist_func_dimension_);
        }
        
        //> perform the normal hsp test
        void HSP_Test(float* query_pointer, std::vector<elementID>& neighbors) const;
        inline void HSP_Test(elementID const index, std::vector<elementID>& neighbors) const {
            HSP_Test(internalDataPointer_(index), neighbors);
        }

        // void constructIndex(std::vector<float> radius_vector) {};

        // void HHSP_Test(){};

        //> print the entire representation for the index
        inline void print_representation(elementID const index) {
            float* index_pointer = internalDataPointer_(index);
            printf("%u: {",(unsigned int) index);
            for (int i = 0; i < (int) dimension_; i++) {
                printf("%.6f,",*index_pointer);
                index_pointer++;
            }
            printf("}\n");
            return;
        }

        
};

#endif // HierarchicalHSP_hpp