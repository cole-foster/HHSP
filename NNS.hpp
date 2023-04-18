#ifndef NNS_hpp
#define NNS_hpp
// Cole Foster
// 2022-11-17

#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot-layer.hpp"

namespace NNS {
    void NNS_2L(unsigned int const queryIndex, PivotLayer& pivotLayer, SparseMatrix &sparseMatrix, unsigned int& nearestNeighbor);
    //void NNS_3L(unsigned int const queryIndex, PivotLayer& pivotLayer, SparseMatrix &sparseMatrix, unsigned int& nearestNeighbor);

    void NNS_BruteForce(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix, unsigned int &nearestNeighbor);
}

#endif // NNS_hpp