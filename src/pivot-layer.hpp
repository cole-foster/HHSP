#ifndef PivotLayer_hpp
#define PivotLayer_hpp

#include <vector> 
#include <tsl/sparse_map.h>
#include <tsl/sparse_set.h>

struct PivotLayer {
   public:
    PivotLayer(){};
    PivotLayer(int const pivotLayerID, int const numberOfPivotLayers, float const radius)
        : pivotLayerID(pivotLayerID), numberOfPivotLayers(numberOfPivotLayers), radius(radius){};
    ~PivotLayer(){};

    //> Add a pivot to this layer
    inline void addPivot(unsigned int pivotIndex) {
        (*pivotIndices).insert(pivotIndex);

        if (pivotLayerID < numberOfPivotLayers - 2) {  // just not bottom layer
            (*pivotChildren)[pivotIndex] = std::vector<unsigned int>{};
            (*maxChildDistance)[pivotIndex] = 0.0f;
        }
    }

    //> Add a pivot to this layer
    inline void addParent(unsigned int pivotIndex, unsigned int parentIndex) {
        (*parentIndices)[pivotIndex] = parentIndex;
    }

    //> Add a child to the pivot in the layer
    inline void addChild(unsigned int const pivotIndex, unsigned int const childIndex) {
        (*pivotChildren)[pivotIndex].push_back(childIndex);
    }

    //> Given this distance to a child/grandchild, update the max distance
    inline void updateMaxChildDistance(unsigned int const pivotIndex, float const distance) {
        if (distance > (*maxChildDistance)[pivotIndex]) {
            (*maxChildDistance)[pivotIndex] = distance;
        }
    }

    //> Layer Information
    int pivotLayerID = 0;
    int numberOfPivotLayers = 1;
    float radius = 0.0f;

    //> Pivot Information
    inline std::shared_ptr<tsl::sparse_set<unsigned int>> const& get_pivotIndices_ptr() const { return pivotIndices; }
    inline std::vector<unsigned int> const& get_pivotChildren(unsigned int const pivotIndex) const {
        return (*pivotChildren)[pivotIndex];
    }
    inline unsigned int const& get_parentIndex(unsigned int const pivotIndex) const {
        return (*parentIndices)[pivotIndex];
    }
    inline float const& get_maxChildDistance(unsigned int const pivotIndex) const {
        return (*maxChildDistance)[pivotIndex];
    }

    // printing/debuggins
    inline unsigned int get_num_pivots() const {
        return (unsigned int) pivotIndices->size();
    } 
    inline float get_radius() const {
        return radius;
    } 


    private:

    // list of all pivots in the layer
    std::shared_ptr<tsl::sparse_set<unsigned int>> pivotIndices = std::make_shared<tsl::sparse_set<unsigned int>> ();

    // map for each pivot to all children
    std::shared_ptr < tsl::sparse_map<unsigned int, std::vector<unsigned int>>> pivotChildren =
        std::make_shared< tsl::sparse_map<unsigned int, std::vector<unsigned int>>> ();

    // map for parent of each pivot
    std::shared_ptr < tsl::sparse_map<unsigned int, unsigned int>> parentIndices =
        std::make_shared<tsl::sparse_map<unsigned int, unsigned int>> ();

    // map for max distance to a child for each pivot
    std::shared_ptr < tsl::sparse_map<unsigned int, float>> maxChildDistance =
        std::make_shared<tsl::sparse_map<unsigned int, float>> ();


};

#endif // PivotLayer_hpp