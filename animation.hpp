#ifndef Animation_hpp
#define Animation_hpp
// Cole Foster
// 2023-04-06

#include <fstream>

struct Animation {
   public:
    Animation(std::string save_directory, int numberOfLayers) : base_path(save_directory), numLayers(numberOfLayers) {
        animation_nns_timestamp.resize(numLayers);
        animation_nns_pivots.resize(numLayers);

        //> HSP
        // HSP search procedure
        animation_hsp_search_iteration.resize(numLayers);
        animation_hsp_search_timestamp.resize(numLayers);
        animation_hsp_search_pivots.resize(numLayers);

        // HSP elimination procedure
        animation_hsp_iteration.resize(numLayers);
        animation_hsp_timestamp.resize(numLayers);
        animation_hsp_pivots.resize(numLayers);
        animation_hsp_value.resize(numLayers);
    };
    ~Animation(){};
    std::string base_path;
    int numLayers;

    void save_pivots(float*& dataPointer, unsigned int const dimension, std::vector<unsigned int>& pivots,
                     int layerID) {
        std::string pivots_out_path =
            std::string(base_path).append("pivots-").append(std::to_string(layerID + 1)).append(".txt");
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
        std::string pivots_out_path = std::string(base_path).append("dataset.txt");
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

    unsigned int timestamp = 0;
    int currentLayer = 0;

    // nns
    std::vector<std::vector<unsigned int>> animation_nns_timestamp{};
    std::vector<std::vector<unsigned int>> animation_nns_pivots{};
    void add_nns_point(unsigned int const index, int layer) {
        for (int i = layer; i < numLayers; i++) {
            animation_nns_timestamp[i].push_back(timestamp);
            animation_nns_pivots[i].push_back(index);
        }
    }
    void save_nns(int layerID) {
        std::string out_path = std::string(base_path).append("nns-").append(std::to_string(layerID + 1)).append(".txt");
        FILE* file = fopen(out_path.c_str(), "w");
        if (file != NULL) {
            fprintf(file, "timestamp,pivotID\n");
            int numpivs = animation_nns_timestamp[layerID].size();
            printf("Layer: %d, num_pivs = %u\n", layerID, numpivs);
            for (int p = 0; p < numpivs; p++) {
                fprintf(file, "%u,%u\n", animation_nns_timestamp[layerID][p], animation_nns_pivots[layerID][p]);
            }
        }
        fclose(file);
    }

    //================================================
    //> HSP ANIMATION
    //================================================
    int iteration = 0;

    //> Storing the HSP Neighbors
    std::vector<unsigned int> animation_hsp_neighbors{};
    void save_neighbors() {
        std::string out_path = std::string(base_path).append("neighbors.txt");
        FILE* file = fopen(out_path.c_str(), "w");
        if (file != NULL) {
            for (int p = 0; p < animation_hsp_neighbors.size(); p++) {
                fprintf(file, "%u\n", animation_hsp_neighbors[p]);
            }
        }
        fclose(file);
    }


    //> Iterations of the HSP Search Procedure
    unsigned int timestamp_search = 0;
    std::vector<std::vector<int>> animation_hsp_search_iteration{};
    std::vector<std::vector<int>> animation_hsp_search_timestamp{};
    std::vector<std::vector<unsigned int>> animation_hsp_search_pivots{};

    /**
     * @brief
     *
     * @param layer      layer of the pivot
     * @param index      id of the pivot
     * @param value      value of outcome: 0 -> safe, 1 -> eliminated, 2 -> intermediate, 3 -> domain examined
     */
    void add_hsp_search_point(int layer, unsigned int index) {
        for (int i = layer; i < numLayers; i++) {
            animation_hsp_search_iteration[i].push_back(iteration);
            animation_hsp_search_timestamp[i].push_back(timestamp_search);
            animation_hsp_search_pivots[i].push_back(index);
        }
    }
    void save_hsp_search(int layerID) {
        std::string out_path =
            std::string(base_path).append("hsp-search-").append(std::to_string(layerID + 1)).append(".txt");
        FILE* file = fopen(out_path.c_str(), "w");
        if (file != NULL) {
            fprintf(file, "iteration,timestamp,pivotID\n");
            int numpivs = animation_hsp_search_iteration[layerID].size();
            printf("Layer: %d, num_pivs = %u\n", layerID, numpivs);
            for (int p = 0; p < numpivs; p++) {
                fprintf(file, "%d,%u,%u\n", animation_hsp_search_iteration[layerID][p],
                        animation_hsp_search_timestamp[layerID][p], animation_hsp_search_pivots[layerID][p]);
            }
        }
        fclose(file);
    }

    // elimination iterations
    unsigned int timestamp_hsp = 0;
    std::vector<std::vector<int>> animation_hsp_iteration{};
    std::vector<std::vector<unsigned int>> animation_hsp_timestamp{};
    std::vector<std::vector<unsigned int>> animation_hsp_pivots{};
    std::vector<std::vector<int>> animation_hsp_value{};

    /**
     * @brief
     *
     * @param layer      layer of the pivot
     * @param index      id of the pivot
     * @param value      value of outcome: 0 -> safe, 1 -> eliminated, 2 -> intermediate, 3 -> domain examined
     */
    void add_hsp_point(int layer, unsigned int index, int value) {
        animation_hsp_iteration[layer].push_back(iteration);
        animation_hsp_timestamp[layer].push_back(timestamp_hsp);
        animation_hsp_pivots[layer].push_back(index);
        animation_hsp_value[layer].push_back(value);
    }
    void save_hsp(int layerID) {
        std::string out_path = std::string(base_path).append("hsp-").append(std::to_string(layerID + 1)).append(".txt");
        FILE* file = fopen(out_path.c_str(), "w");
        if (file != NULL) {
            fprintf(file, "iteration,timestamp,pivotID,outcome\n");
            int numpivs = animation_hsp_timestamp[layerID].size();
            printf("Layer: %d, num_pivs = %u\n", layerID, numpivs);
            for (int p = 0; p < numpivs; p++) {
                fprintf(file, "%d,%u,%u,%d\n", animation_hsp_iteration[layerID][p], animation_hsp_timestamp[layerID][p], 
                        animation_hsp_pivots[layerID][p], animation_hsp_value[layerID][p]);
            }
        }
        fclose(file);
    }
};

#endif  // Animation_hpp