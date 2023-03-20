#ifndef PivotSelection_hpp
#define PivotSelection_hpp

// Cole Foster
// 2022-11-17
// Pivot Selection'
#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot-layer.hpp"

namespace PivotSelection {
    void selectPivots(float radius, SparseMatrix& sparseMatrix, PivotLayer& pivotLayer);

    bool validatePivotSelection(PivotLayer &pivotLayer, SparseMatrix &sparseMatrix);
    void printSet(tsl::sparse_set<unsigned int> const& set);
}

#endif // PivotSelection_hpp