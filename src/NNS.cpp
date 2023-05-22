#include "NNS.hpp"
#include <queue>

namespace NNS {
    //> priority queue
    typedef std::priority_queue<std::pair<float, unsigned int>, std::vector<std::pair<float, unsigned int>>,
                            std::greater<std::pair<float, unsigned int>>> PriorityQueue;

    //> Helper Functions
    float const computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);
    float const getDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);
    void printSet(std::vector<unsigned int> const &set);
};

/**
 * ==============================================================================
 *
 *                                  3-Layer
 *
 * ==============================================================================
 */


/**
 * @brief Get NN by pivot index
 *
 * @param queryIndex
 * @param pivotLayers
 * @param sparseMatrix
 * @param nearestNeighbor
 */
void NNS::NNS_3L(unsigned int const queryIndex, std::vector<PivotLayer> &pivotLayers, SparseMatrix &sparseMatrix,
                    unsigned int &nearestNeighbor) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    sparseMatrix._clear(); // remove all stored distances
    sparseMatrix._addNewReference(queryIndex); // add vector for query distance storage
    tsl::sparse_set<unsigned int> const &pivotIndices = *pivotLayers[0].get_pivotIndices_ptr();
    float dmin = HUGE_VAL;

    PriorityQueue Q1, Q2;
    std::pair<float, unsigned int> topOfQueue;

    //> 1. Iterate Through Top Layer Pivots To Constrain Search Radius, Create Priority Queue
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex1 = (*it1);
        float const distance_Q1 = getDistance(queryIndex, pivotIndex1, sparseMatrix);
        float const radius1 = pivotLayers[0].get_maxChildDistance(pivotIndex1); 

        // get an estimate of the NN by only using pivots, defines search radius
        if (distance_Q1 < dmin) {
            dmin = distance_Q1;
            nearestNeighbor = pivotIndex1;
        }

        // add to queue if it may contain a closer point
        if (distance_Q1 <= dmin + radius1) {
            Q1.push(std::make_pair(distance_Q1,pivotIndex1));
        }
    }

    //> 2. Search Domains of Top Layer Pivots
    while (Q1.size() > 0) {
        topOfQueue = Q1.top();
        Q1.pop();
        unsigned int const pivotIndex1 = topOfQueue.second;
        float const distance_Q1 = topOfQueue.first;
        float const radius1 = pivotLayers[0].get_maxChildDistance(pivotIndex1);

        // Can this pivot contain a closer pivot in the second layer?
        if (distance_Q1 <= dmin + radius1) {
            std::vector<unsigned int> const &pivotDomain = pivotLayers[0].get_pivotChildren(pivotIndex1);

            // Iterate through the pivot domain on second layer
            std::vector<unsigned int>::const_iterator it2;
            for (it2 = pivotDomain.begin(); it2 != pivotDomain.end(); it2++) {
                unsigned int const pivotIndex2 = (*it2);
                float const distance_Q2 = getDistance(queryIndex, pivotIndex2, sparseMatrix);
                float const radius2 = pivotLayers[1].get_maxChildDistance(pivotIndex2);

                // Get an estimate of the NN by only using pivots, defines search radius
                if (distance_Q2 < dmin) {
                    dmin = distance_Q2;
                    nearestNeighbor = pivotIndex2;
                }

                // Add the domain if it may contain a closer point
                if (distance_Q2 <= dmin + radius2) {
                    Q2.push(std::make_pair(distance_Q2, pivotIndex2));
                }
            }
        }
    }

    //> 3. Search Domains of Second Layer Pivots
    while (Q2.size() > 0) {
        topOfQueue = Q2.top();
        Q2.pop();
        unsigned int const pivotIndex2 = topOfQueue.second;
        float const distance_Q2 = topOfQueue.first;
        float const radius2 = pivotLayers[1].get_maxChildDistance(pivotIndex2);

        // Can this pivot contain a closer pivot in the second layer?
        if (distance_Q2 <= dmin + radius2) {
            std::vector<unsigned int> const &pivotDomain = pivotLayers[1].get_pivotChildren(pivotIndex2);

            // iterate through the domain
            std::vector<unsigned int>::const_iterator it3;
            for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                unsigned int const pivotIndex3 = (*it3);
                float const distance_Q3 = getDistance(queryIndex, pivotIndex3, sparseMatrix);

                // Update dmin
                if (distance_Q3 < dmin) {
                    dmin = distance_Q3;
                    nearestNeighbor = pivotIndex3;
                }
            }
        }
    }

    return;
}



/**
 * ==============================================================================
 *
 *                                  2-Layer
 *
 * ==============================================================================
 */

/**
 * @brief Get NN by a 2-Layer Pivot Index
 *
 * @param queryIndex
 * @param pivotLayers
 * @param sparseMatrix
 * @param nearestNeighbor
 */
void NNS::NNS_2L(unsigned int const queryIndex, std::vector<PivotLayer>& pivotLayers, SparseMatrix &sparseMatrix,
                    unsigned int &nearestNeighbor) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    sparseMatrix._clear(); // remove all stored distances
    sparseMatrix._addNewReference(queryIndex); // add vector for query distance storage
    tsl::sparse_set<unsigned int> const &pivotIndices = *pivotLayers[0].get_pivotIndices_ptr();

    // define search radius as distance to closest pivot
    float searchRadius = HUGE_VAL;
    PriorityQueue topPivotsQueue;
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        float const distance_Qp = getDistance(queryIndex, pivotIndex, sparseMatrix);
        topPivotsQueue.push(std::make_pair(distance_Qp,pivotIndex));

        // get an estimate of the NN by only using pivots, defines search radius
        if (distance_Qp < searchRadius) {
            searchRadius = distance_Qp;
            nearestNeighbor = pivotIndex;
        }
    }
    
    // iterate through the queue
    std::vector<unsigned int>::const_iterator it2;
    std::pair<float, unsigned int> topOfQueue;
    while (!topPivotsQueue.empty()) {
        topOfQueue = topPivotsQueue.top();
        topPivotsQueue.pop();
        unsigned int const pivotIndex = topOfQueue.second;
        float const distance_Qp = topOfQueue.first;
        float const radius = pivotLayers[0].get_maxChildDistance(pivotIndex);

        // search the pivot domain
        if (distance_Qp <= searchRadius + radius) {
            std::vector<unsigned int> const &pivotDomain = pivotLayers[0].get_pivotChildren(pivotIndex);
            for (it2 = pivotDomain.begin(); it2 != pivotDomain.end(); it2++) {
                unsigned int const childIndex = (*it2);
                float const distance_Qc = getDistance(queryIndex, childIndex, sparseMatrix);

                // the new nearest neighbor updates the search radius!
                if (distance_Qc < searchRadius) {
                    searchRadius = distance_Qc;
                    nearestNeighbor = childIndex;
                }
            }
        }
    }

    return;
}










/**
 * ==============================================================================
 *
 *                                  BRUTE FORCE
 *
 * ==============================================================================
 */

/**
 * @brief Find NN by brute force
 *
 * @param queryIndex
 * @param datasetSize
 * @param sparseMatrix
 * @param nearestNeighbor
 */
void NNS::NNS_BruteForce(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
                    unsigned int &nearestNeighbor) {
    float nearestDistance = HUGE_VAL;

    for (unsigned int xi = 0; xi < datasetSize; xi++) {
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
 *                                  HELPER FUNCTIONS
 *
 * ==============================================================================
 */

// compute distance using sparse matrix
float const NNS::computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix) {
    return sparseMatrix._computeDistance(index1, index2);
}

// get or retrieve query distnace
float const NNS::getDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix) {
    return sparseMatrix._getDistance(index1, index2);
}

void NNS::printSet(std::vector<unsigned int> const &set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}