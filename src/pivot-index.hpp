#ifndef PivotIndex_hpp
#define PivotIndex_hpp

// Cole Foster
// 2022-11-17
// Pivot Selection'
#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot-layer.hpp"

namespace PivotIndex {
    void Greedy_2L(float radius, SparseMatrix& sparseMatrix, std::vector<PivotLayer>& pivotLayers);
    void Greedy_3L(std::vector<float> radiusVector, SparseMatrix& sparseMatrix, std::vector<PivotLayer>& pivotLayers);

    bool validatePivotSelection(PivotLayer &pivotLayer, SparseMatrix &sparseMatrix);
    bool validatePivotSelection_3L(std::vector<PivotLayer>& pivotLayers, SparseMatrix& sparseMatrix);
    void printSet(std::vector<unsigned int> const& set);
}

#endif // PivotIndex_hpp