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
    void Greedy_MultiLayer(std::vector<float> radiusVector, SparseMatrix& sparseMatrix, std::vector<Pivot>& pivot_list, std::vector<unsigned int>& pivotsPerLayer);

    // bool validatePivotSelection(std::vector<Pivot>& pivotsList, SparseMatrix &sparseMatrix);
    bool validatePivotSelection(std::vector<Pivot>& pivotsList, int const numberOfLayers, SparseMatrix& sparseMatrix);


    void printSet(std::vector<unsigned int> const &set);
}

#endif // PivotIndex_hpp