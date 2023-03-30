#include "NNS.hpp"

#include <bits/stdc++.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace NNS {

typedef std::priority_queue<std::pair<float, const Pivot*>, std::vector<std::pair<float, const Pivot*>>,
                        std::greater<std::pair<float, const Pivot*>>> PriorityQueue;

void recursiveSearch(PriorityQueue& pqueue, unsigned int const queryIndex, unsigned int &index1,
                     float &dmin, SparseMatrix &sparseMatrix, std::vector<float> &queryDistances);

// void depthFirstSearch(const Pivot *pivot, unsigned int const queryIndex, unsigned int &index1, float &dmin,
//                       SparseMatrix &sparseMatrix, std::vector<float> &queryDistances);

}  // namespace NNS


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

    // Define the search radius and create priority queue
    PriorityQueue pqueue{}; // GOTTA USE A FASTER LIBRARY
    std::vector<Pivot>::const_iterator it1, it2;
    for (it1 = pivotsList.begin(); it1 != pivotsList.end(); it1++) {
        const Pivot *pivot = &(*it1);
        float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);
        pqueue.push(std::make_pair(distance,pivot));

        // get an estimate of the NN by only using pivots, defines search radius
        if (distance < searchRadius) {
            searchRadius = distance;
            nearestNeighbor = pivot->_index;
        }
    }

    // go down, recursive baby
    recursiveSearch(pqueue, queryIndex, nearestNeighbor, searchRadius, sparseMatrix, queryDistances);

    return;
}


/**
 * @brief Perform the recursive search
 * 
 * @param pqueue 
 * @param queryIndex 
 * @param index1 
 * @param dmin 
 * @param sparseMatrix 
 * @param queryDistances 
 */
void NNS::recursiveSearch(PriorityQueue& pqueue, unsigned int const queryIndex, unsigned int &index1,
                          float &dmin, SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {
    PriorityQueue nextQueue{}; 

    // iterate through the remainder of the queue   
    std::pair<float, const Pivot*> topOfQueue;
    std::vector<Pivot>::const_iterator it2;
    while (!pqueue.empty()) {
        topOfQueue = pqueue.top();
        const Pivot* pivot = topOfQueue.second;
        float distance1 = topOfQueue.first;

        // add entire domain if it updates
        if (distance1 <= dmin + pivot->_radius) {
            if (pivot->_childCount > 0) {
                for (it2 = pivot->_pivotDomain.begin(); it2 != pivot->_pivotDomain.end(); it2++) {
                    const Pivot* childPivot = &(*it2);
                    float const distance2 = getQueryDistance(queryIndex, childPivot->_index, sparseMatrix, queryDistances);

                    // update search radius
                    if (distance2 < dmin) {
                        dmin = distance2;
                        index1 = childPivot->_index;
                    }

                    // if there are lower children
                    if (childPivot->_childCount > 1) { // 1, only a parent of itself
                        nextQueue.push(std::make_pair(distance2,childPivot));
                    }
                }
            }
        }

        pqueue.pop();
    }

    // go down to next layer
    if (nextQueue.size() > 0) {
        recursiveSearch(nextQueue, queryIndex, index1, dmin, sparseMatrix, queryDistances);
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


// void NNS::depthFirstSearch(const Pivot *pivot, unsigned int const queryIndex, unsigned int &index1, float &dmin,
//                            SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {
//     // store the closest point
//     const Pivot *closestPivot;
//     float closestPivotDistance = HUGE_VAL;

//     // iterate through domain to find closest point
//     std::vector<Pivot>::const_iterator it1;
//     for (it1 = pivot->_pivotDomain.begin(); it1 != pivot->_pivotDomain.end(); it1++) {
//         const Pivot *childPivot = &(*it1);
//         float const distance = getQueryDistance(queryIndex, childPivot->_index, sparseMatrix, queryDistances);

//         // find the closest child
//         if (distance < closestPivotDistance) {
//             closestPivotDistance = distance;
//             closestPivot = childPivot;
//         }

//         // update dmin
//         if (distance < dmin) {
//             dmin = distance;
//             index1 = childPivot->_index;
//         }
//     }

//     // go down
//     if (closestPivot->_childCount > 0) {
//         depthFirstSearch(closestPivot, queryIndex, index1, dmin, sparseMatrix, queryDistances);
//     }

//     return;
// }