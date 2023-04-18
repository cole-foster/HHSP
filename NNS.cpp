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
 * @brief Get NN by pivot index
 *
 * @param queryIndex
 * @param pivotLayer
 * @param sparseMatrix
 * @param nearestNeighbor
 */
void NNS::NNS_2L(unsigned int const queryIndex, PivotLayer &pivotLayer, SparseMatrix &sparseMatrix,
                    unsigned int &nearestNeighbor) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    sparseMatrix._clear(); // remove all stored distances
    sparseMatrix._addNewReference(queryIndex); // add vector for query distance storage
    tsl::sparse_set<unsigned int> const &pivotIndices = *pivotLayer.get_pivotIndices_ptr();

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
        float const radius = pivotLayer.get_maxChildDistance(pivotIndex);

        // search the pivot domain
        if (distance_Qp <= searchRadius + radius) {
            std::vector<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(pivotIndex);
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