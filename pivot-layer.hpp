#ifndef PivotLayer_hpp
#define PivotLayer_hpp

// Cole Foster
// 2022-11-17
// Pivot Layer for PLH
#include <tsl/sparse_map.h>
#include <tsl/sparse_set.h>

/**
 * @brief Struct used during pivot selection to hold pivot indices in each layer and pivot domain members
 *
 */
struct PivotLayer {
   public:
    PivotLayer(){};
    PivotLayer(int const pivotLayerID, int const numberOfPivotLayers, float const radius)
        : pivotLayerID(pivotLayerID), numberOfPivotLayers(numberOfPivotLayers), radius(radius){};
    ~PivotLayer(){};

    inline void addPivot(unsigned int pivotIndex) {
        (*pivotIndices).insert(pivotIndex);

        if (pivotLayerID < numberOfPivotLayers - 2) {  // just not bottom layer
            (*pivotChildren)[pivotIndex] = tsl::sparse_set<unsigned int>{};
        }
    }
    inline void addPivotChild(unsigned int pivotIndex, unsigned int childIndex) {
        (*pivotChildren)[pivotIndex].insert(childIndex);
    }

    int pivotLayerID = 0;
    int numberOfPivotLayers = 1;
    float radius = 0.0f;

    std::shared_ptr<tsl::sparse_set<unsigned int>> pivotIndices = std::make_shared<tsl::sparse_set<unsigned int>>();
    std::shared_ptr<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>> pivotChildren =
        std::make_shared<tsl::sparse_map<unsigned int, tsl::sparse_set<unsigned int>>>();

    inline std::shared_ptr<tsl::sparse_set<unsigned int>> const& get_pivotIndices_ptr() const { return pivotIndices; }
    inline tsl::sparse_set<unsigned int> const& get_pivotChildren(unsigned int const pivotIndex) const {
        return (*pivotChildren)[pivotIndex];
    }
};

#endif  // PivotLayer_hpp