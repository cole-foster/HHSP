#ifndef GHSP_hpp
#define GHSP_hpp
// Cole Foster
// 2022-11-17

#include <vector> 
#include "sparse-matrix.hpp"
#include "pivot-layer.hpp"

namespace GHSP {

    void Hierarchical_HSP_Test(unsigned int const queryIndex, std::vector<PivotLayer> &coverTree, SparseMatrix &sparseMatrix, std::vector<unsigned int> &neighbors);

    void GHSP_2L(unsigned int const queryIndex, std::vector<PivotLayer> &pivotLayers, SparseMatrix &sparseMatrix, std::vector<unsigned int> &neighbors);
    void GHSP_3L(unsigned int const queryIndex, std::vector<PivotLayer> &pivotLayers, SparseMatrix &sparseMatrix, std::vector<unsigned int> &neighbors);
    
    void HSP(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix, std::vector<unsigned int>& neighbors);
}

#endif // GHSP_hpp