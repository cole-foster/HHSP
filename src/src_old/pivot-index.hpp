#ifndef PivotIndex_hpp
#define PivotIndex_hpp

// Cole Foster
// 2022-11-17
// Pivot Selection'
#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot-layer.hpp"

namespace PivotIndex {
    void CoverTree_Greedy(std::vector<float> radiusVector, SparseMatrix& sparseMatrix, std::vector<PivotLayer>& coverTree);
    bool validateCoverTree(std::vector<PivotLayer>& pivotLayers, SparseMatrix& sparseMatrix);
    void printSet(std::vector<unsigned int> const& set);
}

#endif // PivotIndex_hpp