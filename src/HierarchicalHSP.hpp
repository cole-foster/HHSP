#ifndef HierarchicalHSP_hpp
#define HierarchicalHSP_hpp

#include "metrics.h"

/*
    Heavily influenced by https://github.com/nmslib/hnswlib
*/

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
        char* data_pointer_ = nullptr;  // char for indexing by bytes
        size_t size_data_per_element_ = 0;
        
        // distance metric: from metrics.h, needs void*
        DISTFUNC<float> distance_function_;
        void *dist_func_dimension_ = nullptr;

    public:
        HierarchicalHSP() {};
        HierarchicalHSP(size_t const dimension, size_t const dataset_size, float*& data_pointer): dimension_(dimension),
                dataset_size_(dataset_size) {
            distance_function_ = metrics::EuclideanDistance;
            dist_func_dimension_ = &dimension_;
            data_pointer_ = reinterpret_cast<char*> (data_pointer); // char: 1 byte per element, float: 4 bytes per el
        };
        ~HierarchicalHSP(){};

        // void HSP_Test(){};

        // void constructIndex(std::vector<float> radius_vector) {};

        // void HHSP_Test(){};
};

#endif // HierarchicalHSP_hpp