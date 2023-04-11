#ifndef Animation_hpp
#define Animation_hpp
// Cole Foster
// 2023-04-06

#include "datasets.hpp"

class Animation {
    std::string base_path;

    void save_pivots(float*& dataPointer, unsigned int const dimension, std::vector<unsigned int>& pivots) {
        std::string pivots_out_path = std::string(resultsDirectory).append("pivots.txt");
        FILE* file = fopen(pivots_out_path.c_str(), "w");
        if (file != NULL) {
            for (int i = 0; i < (int)pivots.size(); i++) {
                unsigned int pivotIndex = pivots[i];
                fprintf(file, "%u,%.6f,%.6f\n", pivotIndex, dataPointer[pivotIndex * dimension],
                        dataPointer[pivotIndex * dimension + 1]);
            }
        }
        fclose(file);
    }

    void save_dataset(float*& dataPointer, unsigned int const dimension, unsigned int const datasetSize) {
        std::string pivots_out_path = std::string(resultsDirectory).append("dataset.txt");
        FILE* file = fopen(pivots_out_path.c_str(), "w");
        if (file != NULL) {
            for (int i = 0; i < datasetSize; i++) {
                fprintf(file, "%u,%.6f,%.6f\n", i, dataPointer[i * dimension], dataPointer[i * dimension + 1]);
            }
        }
        fclose(file);
    }

    void save_query(float*& dataPointer, unsigned int const dimension, unsigned int const queryID) {
        std::string query_out_path = std::string(base_path).append("query.txt");
        FILE* file = fopen(query_out_path.c_str(), "w");
        if (file != NULL) {
            fprintf(file, "%.6f,%.6f\n", dataPointer[queryID * dimension + 0], dataPointer[queryID * dimension + 1]);
        }
        fclose(file);
    }
}

#endif  // Animation_hpp