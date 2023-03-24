#ifndef GHSP_hpp
#define GHSP_hpp
// Cole Foster
// 2023-03-24

#include <vector>

#include "sparse-matrix.hpp"
#include "pivot.hpp"

namespace GHSP {

void GHSP_Search(unsigned int const queryIndex, std::vector<Pivot> const& pivotsList, SparseMatrix &sparseMatrix,
                       std::vector<unsigned int> &neighbors);
void HSP_Search(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
         std::vector<unsigned int> &neighbors);

// helper functions
float const computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);
float const getQueryDistance(unsigned int const queryIndex, unsigned int const index2, SparseMatrix &sparseMatrix,
                             std::vector<float> &distanceList);
void printSet(std::vector<unsigned int> const &set);
}  // namespace GHSP

#endif  // GHSP_hpp