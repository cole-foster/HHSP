#include "NNS.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace NNS {
void recursiveSearch(std::vector<const Pivot*> pivotList,  unsigned int const queryIndex, unsigned int &index1, float &dmin,
                     SparseMatrix &sparseMatrix, std::vector<float> &queryDistances);
void recursiveSearch1(Pivot const &parentPivot, unsigned int const queryIndex, unsigned int &index1, float &dmin,
                          SparseMatrix &sparseMatrix, std::vector<float> &queryDistances);
}

/**
 * @brief Get NN by pivot index
 *
 * @param queryIndex
 * @param pivotLayer
 * @param sparseMatrix
 * @param nearestNeighbor
 */
void NNS::Search(unsigned int const queryIndex, std::vector<Pivot> const &pivotsList, SparseMatrix &sparseMatrix,
                 unsigned int &nearestNeighbor) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;

    // initializations for query
    float searchRadius = HUGE_VAL;
    std::vector<float> queryDistances{};
    queryDistances.resize(datasetSize, -1.0f);

    // define search radius as distance to closest pivot
    std::vector<Pivot>::const_iterator it1, it2;
    for (it1 = pivotsList.begin(); it1 != pivotsList.end(); it1++) {
        Pivot const &pivot = (*it1);
        float const distance_Qp = getQueryDistance(queryIndex, pivot._index, sparseMatrix, queryDistances);

        // get an estimate of the NN by only using pivots, defines search radius
        if (distance_Qp < searchRadius) {
            searchRadius = distance_Qp;
            nearestNeighbor = pivot._index;
        }
    }

    // find the closest point from those pivot domains within the search radius
    std::vector<const Pivot*> nextPivotList{};
    for (it1 = pivotsList.begin(); it1 != pivotsList.end(); it1++) {
        Pivot const &pivot = (*it1);
        float const distance = getQueryDistance(queryIndex, pivot._index, sparseMatrix, queryDistances);

        // add entire domain if it updates
        if (distance <= searchRadius + pivot._radius) {
            if (pivot._childCount > 0) {
                for (it2 = pivot._pivotDomain.begin(); it2 != pivot._pivotDomain.end(); it2++) {
                    nextPivotList.push_back(&(*it2));
                }
            }
        }
    }

    recursiveSearch(nextPivotList, queryIndex, nearestNeighbor, searchRadius, sparseMatrix, queryDistances);

    return;
}



void NNS::recursiveSearch(std::vector<const Pivot*> pivotList,  unsigned int const queryIndex, unsigned int &index1, float &dmin,
                     SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {

    // update dmin with all pivots in the list
    std::vector<const Pivot*>::const_iterator it1;
    for (it1 = pivotList.begin(); it1 != pivotList.end(); it1++) {
        const Pivot* pivot = (*it1);
        float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

        // update active NN only if its not an HSP neighbor already
        if (distance < dmin) {
            dmin = distance;
            index1 = pivot->_index;
        }
    }

    // collect pivot domains
    std::vector<const Pivot*> nextPivotList{};
    std::vector<Pivot>::const_iterator it2;
    for (it1 = pivotList.begin(); it1 != pivotList.end(); it1++) {
        const Pivot* pivot = (*it1);
        float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

        // add entire domain if it updates
        if (distance <= dmin + pivot->_radius) {
            if (pivot->_childCount > 0) {
                for (it2 = pivot->_pivotDomain.begin(); it2 != pivot->_pivotDomain.end(); it2++) {
                    nextPivotList.push_back(&(*it2));
                }
            }
        }
    }

    // go down to next layer
    if (nextPivotList.size() > 0) {
        recursiveSearch(nextPivotList, queryIndex, index1, dmin, sparseMatrix, queryDistances);
    }

    return;
}


/**
 * @brief Find NN by brute force
 *
 * @param queryIndex
 * @param datasetSize
 * @param sparseMatrix
 * @param nearestNeighbor
 */
void NNS::Search_BF(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
                    unsigned int &nearestNeighbor) {
    float nearestDistance = HUGE_VAL;

    for (unsigned int xi = 0; xi < datasetSize; xi++) {
        if (queryIndex == xi) continue;
        float const distance_Qi = sparseMatrix._computeDistance(queryIndex, xi);
        if (distance_Qi < nearestDistance) {
            nearestDistance = distance_Qi;
            nearestNeighbor = xi;
        }
    }
    return;
}

/**
 * ==============================================================================
 *
 *                      HELPER FUNCTIONS
 *
 * ==============================================================================
 */

void NNS::printSet(std::vector<unsigned int> const &set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}

// compute distance using sparse matrix
float const NNS::computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix) {
    return sparseMatrix._computeDistance(index1, index2);
}

// get or retrieve query distnace
float const NNS::getQueryDistance(unsigned int const queryIndex, unsigned int const index2, SparseMatrix &sparseMatrix,
                                  std::vector<float> &distanceList) {
    float distance = distanceList[index2];
    if (distance < 0) {
        distance = sparseMatrix._computeDistance(queryIndex, index2);
        distanceList[index2] = distance;
    }
    return distance;
}


// /**
//  * @brief Get NN by pivot index
//  *
//  * @param queryIndex
//  * @param pivotLayer
//  * @param sparseMatrix
//  * @param nearestNeighbor
//  */
// void NNS::Search1(unsigned int const queryIndex, std::vector<Pivot> const &pivotsList, SparseMatrix &sparseMatrix,
//                  unsigned int &nearestNeighbor) {
//     unsigned int const datasetSize = sparseMatrix._datasetSize;

//     // initializations for query
//     float searchRadius = HUGE_VAL;
//     std::vector<float> queryDistances{};
//     queryDistances.resize(datasetSize, -1.0f);

//     // define search radius as distance to closest pivot
//     std::vector<Pivot>::const_iterator it1, it2;
//     for (it1 = pivotsList.begin(); it1 != pivotsList.end(); it1++) {
//         Pivot const &pivot = (*it1);
//         float const distance_Qp = getQueryDistance(queryIndex, pivot._index, sparseMatrix, queryDistances);

//         // get an estimate of the NN by only using pivots, defines search radius
//         if (distance_Qp < searchRadius) {
//             searchRadius = distance_Qp;
//             nearestNeighbor = pivot._index;
//         }
//     }

//     // find the closest point from those pivot domains within the search radius
//     for (it1 = pivotsList.begin(); it1 != pivotsList.end(); it1++) {
//         Pivot const &pivot = (*it1);
//         recursiveSearch1(pivot, queryIndex, nearestNeighbor, searchRadius, sparseMatrix, queryDistances);
//     }

//     return;
// }


// /**
//  * @brief
//  *
//  * @param parentPivot
//  * @param queryIndex
//  * @param neighbors
//  * @param index1
//  * @param dmin
//  * @param sparseMatrix
//  * @param queryDistances
//  */
// void NNS::recursiveSearch1(Pivot const &parentPivot, unsigned int const queryIndex, unsigned int &index1, float &dmin,
//                           SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {
//     if (parentPivot._childCount <= 0) return;

//     // iterate through pivot domain
//     std::vector<Pivot>::const_iterator it1;
//     for (it1 = parentPivot._pivotDomain.begin(); it1 != parentPivot._pivotDomain.end(); it1++) {
//         Pivot const &childPivot = (*it1);
//         float const distance_Qc = getQueryDistance(queryIndex, childPivot._index, sparseMatrix, queryDistances);

//         // update active NN only if its not an HSP neighbor already
//         if (distance_Qc < dmin) {
//             dmin = distance_Qc;
//             index1 = childPivot._index;
//         }

//         // check if domain can hold a closer point
//         if (distance_Qc <= dmin + childPivot._maxChildDistance) {
//             if (childPivot._childCount > 0) {
//                 recursiveSearch1(childPivot, queryIndex, index1, dmin, sparseMatrix, queryDistances);
//             }
//         }
//     }

//     return;
// }