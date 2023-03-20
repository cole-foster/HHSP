#ifndef GHSP_hpp
#define GHSP_hpp
// Cole Foster
// 2022-11-17

#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot-layer.hpp"

namespace GHSP {

    // NNS
    void pivotNNS(unsigned int const queryIndex, PivotLayer& pivotLayer, SparseMatrix &sparseMatrix, unsigned int& nearestNeighbor);
    void GHSP_Search(unsigned int const queryIndex, PivotLayer &pivotLayer, SparseMatrix &sparseMatrix, std::vector<unsigned int> &neighbors);
    
    // brute-force approaches
    void bruteNNS(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix, unsigned int &nearestNeighbor);
    void HSP(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix, std::vector<unsigned int>& neighbors);

    // helper functions
    float const computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix& sparseMatrix);
    float const getQueryDistance(unsigned int const queryIndex, unsigned int const index2, SparseMatrix& sparseMatrix, std::vector<float>& distanceList);
    void printSet(std::vector<unsigned int> const &set);
}

#endif // GHSP_hpp