#ifndef NNS_hpp
#define NNS_hpp
// Cole Foster
// 2023-03-24

#include <vector>

#include "pivot.hpp"
#include "sparse-matrix.hpp"
#include "animation.hpp"

// Class for performing NNS on a pivot-index
namespace NNS {
void Search_BF(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
               unsigned int &nearestNeighbor);
void Search(unsigned int const queryIndex, std::vector<Pivot> const &pivotsList, SparseMatrix &sparseMatrix,
            unsigned int &nearestNeighbor, Animation& anim);
// void Search1(unsigned int const queryIndex, std::vector<Pivot> const &pivotsList, SparseMatrix &sparseMatrix,
//             unsigned int &nearestNeighbor);

// helper functions
float const computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);
float const getDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);
float const getQueryDistance(unsigned int const queryIndex, unsigned int const index2, SparseMatrix &sparseMatrix,
                             std::vector<float> &distanceList);
void printSet(std::vector<unsigned int> const &set);
}  // namespace NNS

#endif  // NNS_hpp