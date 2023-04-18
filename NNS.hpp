#ifndef NNS_hpp
#define NNS_hpp
// Cole Foster
// 2022-11-17

#include <vector>

#include "pivot-layer.hpp"
#include "sparse-matrix.hpp"

namespace NNS {
void NNS_2L(unsigned int const queryIndex, std::vector<PivotLayer>& pivotLayers, SparseMatrix& sparseMatrix,
            unsigned int& nearestNeighbor);
void NNS_3L(unsigned int const queryIndex, std::vector<PivotLayer>& pivotLayers, SparseMatrix& sparseMatrix,
            unsigned int& nearestNeighbor);

void NNS_BruteForce(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix& sparseMatrix,
                    unsigned int& nearestNeighbor);
}  // namespace NNS

#endif  // NNS_hpp