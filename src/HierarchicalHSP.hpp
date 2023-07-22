#ifndef HierarchicalHSP_hpp
#define HierarchicalHSP_hpp

#include <cstdio>
#include <vector>
#include <algorithm>
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
        float* data_pointer_ = nullptr;

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
        void HSP_Test(float* query_pointer, std::vector<elementID>& neighbors) const {
            neighbors.clear();

            // only perform on k closest elements
            std::vector<std::pair<float, elementID>> active_list{};
            active_list.reserve(dataset_size_);

            // find next nearest neighbor and create list of distances
            elementID index1;
            float distance_Q1 = HUGE_VAL;
            for (elementID index = 0; index < (elementID) dataset_size_; index++) {
                float const distance = compute_distance(query_pointer, internalDataPointer_(index));
                if (distance <= 0) {
                    neighbors.push_back(index);
                    continue;
                }
                if (distance < distance_Q1) {
                    distance_Q1 = distance;
                    index1 = index;
                }
                active_list.emplace_back(distance, index);
            }

            // perform the hsp loop
            while (active_list.size() > 0) {
                neighbors.push_back(index1);
                float* index1_ptr = internalDataPointer_(index1);

                //  - set up for the next hsp neighbor
                elementID index1_next;
                float distance_Q1_next = HUGE_VAL;

                //  - initialize the active_list for next iteration
                std::vector<std::pair<float, elementID>> active_list_copy = active_list;
                active_list.clear();

                //  - check each point for elimination
                for (int it2 = 0; it2 < (int)active_list_copy.size(); it2++) {
                    elementID const index2 = active_list_copy[it2].second;
                    float const distance_Q2 = active_list_copy[it2].first;
                    if (index2 == index1) continue;
                    float const distance_12 = compute_distance(index1_ptr, internalDataPointer_(index2));

                    // check the hsp inequalities: add if not satisfied
                    if (distance_Q1 >= distance_Q2 || distance_12 >= distance_Q2) {
                        active_list.emplace_back(distance_Q2, index2);
                        if (distance_Q2 < distance_Q1_next) {
                            distance_Q1_next = distance_Q2;
                            index1_next = index2;
                        }
                    }
                }

                // setup the next hsp neighbor
                index1 = index1_next;
                distance_Q1 = distance_Q1_next;
            }

            return;
        }
        inline void HSP_Test(elementID const index, std::vector<elementID>& neighbors) const {
            HSP_Test(internalDataPointer_(index), neighbors);
        }

        // void constructIndex(std::vector<float> radius_vector) {};

        // void HHSP_Test(){};

        //> print the entire representation for the index
        inline void print_representation(elementID const index) {
            float* index_pointer = internalDataPointer_(index);
            printf("%u: [",(unsigned int) index);
            for (int i = 0; i < (int) dimension_; i++) {
                printf("%.6f,",*index_pointer);
                index_pointer++;
            }
            printf("]\n");
            return;
        }

        inline void print_vector(std::vector<elementID>& vec) const {
            printf("[");
            for (int i = 0; i < (int) vec.size(); i++) {
                printf("%u, ",vec[i]);
            }
            printf("]\n");
            return;
        }

        
};

#endif // HierarchicalHSP_hpp