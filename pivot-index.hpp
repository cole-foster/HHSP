#ifndef PivotIndex_hpp
#define PivotIndex_hpp
// Cole Foster
// 2023-03-20
// Constructing an index of pivots/pivot domains

#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot.hpp"

namespace PivotIndex {
    void Greedy2L(float radius, SparseMatrix& sparseMatrix, std::vector<Pivot>& pivotsList);

    bool validatePivotSelection(std::vector<Pivot>& pivotsList, SparseMatrix &sparseMatrix);
}

#endif // PivotIndex_hpp