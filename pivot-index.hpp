#ifndef PivotIndex_hpp
#define PivotIndex_hpp

// Cole Foster
// 2022-11-17
// Pivot Selection'
#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot-layer.hpp"

namespace PivotIndex {
    void Greedy_2L(float radius, SparseMatrix& sparseMatrix, PivotLayer& pivotLayer);
    void Greedy_3L(std::vector<float> radiusVector, SparseMatrix& sparseMatrix, std::vector<PivotLayer>& pivotLayer);

    bool validatePivotSelection(PivotLayer &pivotLayer, SparseMatrix &sparseMatrix);
    void printSet(tsl::sparse_set<unsigned int> const& set);
}

#endif // PivotIndex_hpp