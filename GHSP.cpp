#include "GHSP.hpp"

#include <cmath>
#include <cstdio>

// void GHSP::GHSP_Search(unsigned int const queryIndex, PivotLayer &pivotLayer, SparseMatrix &sparseMatrix,
//                        std::vector<unsigned int> &neighbors) {
//     neighbors.clear();

//     // initialize query distance storage
//     unsigned int datasetSize = sparseMatrix._datasetSize;
//     std::vector<float> queryDistances{};
//     queryDistances.resize(datasetSize, -1.0f);

//     // everything that can potentially become an HSP neighbor
//     tsl::sparse_set<unsigned int> const &pivotIndices = *pivotLayer.get_pivotIndices_ptr();
//     float const radius = pivotLayer.radius;
//     std::vector<unsigned int> A1(pivotIndices.begin(), pivotIndices.end());  // active top layer pivots
//     std::vector<unsigned int> I1{};                                          // intermediate top layer pivots
//     std::vector<unsigned int> A2{};                                          // active bottom layer points
//     std::vector<unsigned int> I2{};                                          // intermeidate bottom layer points

//     // find the HSP neighbors in a loop
//     while (A1.size() > 0 || I1.size() > 0 || A2.size() > 0) {

//         /**
//          * ============================================================
//          * @brief STEP 1: Collect dmin as preliminary closest active point
//          * ============================================================
//          */
//         unsigned int index1;  // next hsp neighbor
//         float dmin = HUGE_VAL;

//         // update dmin with closest active point in A2
//         for (unsigned int it2 = 0; it2 < A2.size(); it2++) {
//             unsigned int const x2 = A2[it2];
//             float const distance_Q2 = getQueryDistance(queryIndex, x2, sparseMatrix, queryDistances);
//             if (distance_Q2 < dmin) {
//                 dmin = distance_Q2;
//                 index1 = x2;
//             }
//         }

//         // update dmin with closest active pivot in A1
//         for (unsigned int it2 = 0; it2 < A1.size(); it2++) {
//             unsigned int const p2 = A1[it2];
//             float const distance_Q2 = getQueryDistance(queryIndex, p2, sparseMatrix, queryDistances);
//             if (distance_Q2 < dmin) {
//                 dmin = distance_Q2;
//                 index1 = p2;
//             }
//         }

//         // update dmin with only validated pivots in I1
//         for (unsigned int it2 = 0; it2 < I1.size(); it2++) {
//             unsigned int const p2 = I1[it2];
//             float const distance_Q2 = getQueryDistance(queryIndex, p2, sparseMatrix, queryDistances);

//             // if p2 an hsp neighbor, ignore
//             if (std::find(neighbors.begin(), neighbors.end(), p2) != neighbors.end()) continue;

//             // can update dmin if validated
//             if (distance_Q2 < dmin) {

//                 // validate against existing hsp neighbors
//                 bool flag_validated = true;
//                 for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
//                     unsigned int const x1_bar = neighbors[it1];
//                     float const distance_Q1 = getQueryDistance(queryIndex, x1_bar, sparseMatrix, queryDistances);
//                     float const distance_12 = computeDistance(x1_bar, p2, sparseMatrix);
//                     if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
//                         flag_validated = false;
//                         break;
//                     }
//                 }

//                 // update dmin if p2 is valid
//                 if (flag_validated) {
//                     dmin = distance_Q2;
//                     index1 = p2;
//                 }
//             }
//         }

//         /**
//          * ============================================================
//          * @brief STEP 2: Search the domains of pivots that could contain nearest neighbor
//          * ============================================================
//          */

//         // update dmin with domains of active pivots
//         for (unsigned int it2 = 0; it2 < A1.size(); it2++) {
//             unsigned int const p2 = A1[it2];
//             float const distance_Qp2 = getQueryDistance(queryIndex, p2, sparseMatrix, queryDistances);

//             // could pivot domain contain a closer point?
//             if (distance_Qp2 <= dmin + radius) {
//                 tsl::sparse_set<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);

//                 // iterate through the pivot domain
//                 tsl::sparse_set<unsigned int>::const_iterator it3;
//                 for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
//                     unsigned int const x2 = (*it3);
//                     float const distance_Qx2 = getQueryDistance(queryIndex, x2, sparseMatrix, queryDistances);

//                     // this can update dmin!
//                     if (distance_Qx2 < dmin) {
//                         dmin = distance_Qx2;
//                         index1 = x2;
//                     }
//                 }
//             }
//         }

//         // update dmin with domains of intermediate pivots
//         for (std::vector<unsigned int>::iterator it2 = I1.begin(); it2 != I1.end(); /* iterate in loop */) {
//             unsigned int const p2 = (*it2);
//             float const distance_Qp2 = getQueryDistance(queryIndex, p2, sparseMatrix, queryDistances);

//             // could pivot domain contain a closer point?
//             if (distance_Qp2 > dmin + radius) {
//                 ++it2;  // iterate

//             } else {
//                 it2 = I1.erase(it2);  // entire domain needs to be examined. Might as well remove pivot from list

//                 // iterate through domain
//                 tsl::sparse_set<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);
//                 tsl::sparse_set<unsigned int>::const_iterator it3;
//                 for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
//                     unsigned int const x2 = (*it3);
//                     float const distance_Q2 = getQueryDistance(queryIndex, x2, sparseMatrix, queryDistances);

//                     // if child an hsp neighbor, ignore
//                     if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;

//                     // validate point against all previous HSP neighbors
//                     bool flag_validated = true;
//                     for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
//                         unsigned int const x1_bar = neighbors[it1];
//                         float const distance_Q1 = getQueryDistance(queryIndex, x1_bar, sparseMatrix, queryDistances);
//                         float const distance_12 = computeDistance(x1_bar, x2, sparseMatrix);
//                         if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
//                             flag_validated = false;
//                             break;
//                         }
//                     }

//                     // since not invalidated, add to active points list
//                     if (flag_validated) {
//                         A2.push_back(x2);

//                         // this can update dmin!
//                         if (distance_Q2 < dmin) {
//                             dmin = distance_Q2;
//                             index1 = x2;
//                         }
//                     }
//                 }
//             }
//         }
//         if (dmin > 10000) continue; // large, magic number signaling no neighbor possible

//         /**
//          * ============================================================
//          * @brief STEP 3: Perform GHSP test on all pivots agasinst GHSP between Q,x1
//          * ============================================================
//          */

//         unsigned int const x1 = index1;
//         float const distance_Q1 = getQueryDistance(queryIndex, x1, sparseMatrix, queryDistances);
//         neighbors.push_back(x1);

//         // Iterate through intermediate pivots
//         // Three Options:
//         //  1. Entire Pivot Domain Invalidated
//         //  2. Pivot Remains in List
//         //  3. Entire Domain must be Examined by Brute Force
//         I2.clear();
//         for (std::vector<unsigned int>::iterator it2 = I1.begin(); it2 != I1.end(); /* iterate in loop */) {
//             unsigned int const p2 = (*it2);
//             float const distance_Q2 = getQueryDistance(queryIndex, p2, sparseMatrix, queryDistances);

//             // Pivot cannot be excluded by x1 or any future HSP neighbor
//             // Domain must be validated against all prior HSP neighbor
//             if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
//                 tsl::sparse_set<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);
//                 I2.insert(I2.end(), pivotDomain.begin(), pivotDomain.end());
//                 it2 = I1.erase(it2);

//             } else {
//                 float const distance_12 = computeDistance(x1, p2, sparseMatrix);

//                 // Proposition: Entire Domain Invalidated
//                 if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
//                     it2 = I1.erase(it2);
//                 }
//                 // Else, Remains in List
//                 else {
//                     ++it2;
//                 }
//             }
//         }

//         // Iterate through active pivots
//         // Four Options:
//         //  1. Entire Pivot Domain Invalidated
//         //  2. Pivot Remains in List
//         //  3. Pivot Moves to Intermediate Pivot List
//         //  4. Entire Domain must be Examined by Brute Force

//         // iterate through active pivots list
//         for (std::vector<unsigned int>::iterator it2 = A1.begin(); it2 != A1.end(); /* iterate in loop */) {
//             unsigned int const p2 = (*it2);
//             float const distance_Q2 = getQueryDistance(queryIndex, p2, sparseMatrix, queryDistances);

//             // Pivot cannot be excluded by x1 or any future HSP neighbor
//             if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
//                 tsl::sparse_set<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);
//                 A2.insert(A2.end(), pivotDomain.begin(), pivotDomain.end());
//                 it2 = A1.erase(it2);

//             } else {
//                 float const distance_12 = computeDistance(x1, p2, sparseMatrix);

//                 // Proposition: Entire Domain Invalidated
//                 if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
//                     it2 = A1.erase(it2);
//                 }
//                 // Proposition: Entire Domain Safe
//                 else if (distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
//                     ++it2;
//                 }
//                 // Otherwise, Pivot Falls in Intermediate Region
//                 else {
//                     I1.push_back(p2);
//                     it2 = A1.erase(it2);
//                 }
//             }
//         }

//         /**
//          * ============================================================
//          * @brief STEP 4: Perform normal HSP test on all points in with HSP(Q,x1)
//          * ============================================================
//          */
//         for (std::vector<unsigned int>::iterator it2 = A2.begin(); it2 != A2.end(); /* iterate in loop */) {
//             unsigned int const x2 = (*it2);
//             // remove HSP neighbors from this list, never
//             if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) {
//                 it2 = A2.erase(it2);
//                 continue;
//             }
//             float const distance_Q2 = getQueryDistance(queryIndex, x2, sparseMatrix, queryDistances);
//             float const distance_12 = computeDistance(x1, x2, sparseMatrix);

//             // HSP test
//             if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
//                 it2 = A2.erase(it2);
//             } else {
//                 ++it2;
//             }
//         }

//         /**
//          * ============================================================
//          * @brief STEP 7: Validate all points in I2 against all HSP neighbors
//          * ============================================================
//          */
//         for (unsigned int it2 = 0; it2 < I2.size(); it2++) {
//             unsigned int const x2 = I2[it2];

//             // remove HSP neighbors from this list, never
//             if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;

//             float const distance_Q2 = getQueryDistance(queryIndex, x2, sparseMatrix, queryDistances);

//             // validate against existing hsp neighbors
//             bool flag_validated = true;
//             for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
//                 unsigned int const x1_bar = neighbors[it1];
//                 float const distance_Q1 = getQueryDistance(queryIndex, x1_bar, sparseMatrix, queryDistances);
//                 float const distance_12 = computeDistance(x1_bar, x2, sparseMatrix);
//                 if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
//                     flag_validated = false;
//                     break;
//                 }
//             }

//             // add as an active point if validated
//             if (flag_validated) {
//                 A2.push_back(x2);
//             }
//         }
//     }

//     return;
// }

/**
 * @brief Get NN by pivot index
 *
 * @param queryIndex
 * @param pivotLayer
 * @param sparseMatrix
 * @param nearestNeighbor
 */
void GHSP::pivotNNS(unsigned int const queryIndex, std::vector<Pivot> const &pivotsList, SparseMatrix &sparseMatrix,
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
    for (it1 = pivotsList.begin(); it1 != pivotsList.end(); it1++) {
        Pivot const &pivot = (*it1);
        float const radius = pivot._maxChildDistance;
        float const distance_Qp =
            getQueryDistance(queryIndex, pivot._index, sparseMatrix, queryDistances);

        // condition where pivot may contain nn
        if (distance_Qp <= searchRadius + radius) {

            std::vector<Pivot> const& pivotDomain = pivot._pivotDomain;
            for (it2 = pivotDomain.begin(); it2 != pivotDomain.end(); it2++) {
                unsigned int const childIndex = (*it2)._index;
                float const distance_Qc = getQueryDistance(queryIndex, childIndex, sparseMatrix, queryDistances);

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
 *                      BRUTE FORCE
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
void GHSP::bruteNNS(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
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
 * @brief perform normal HSP Test
 *
 * @param queryIndex
 * @param datasetSize
 * @param sparseMatrix
 * @param neighbors
 */
void GHSP::HSP(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
               std::vector<unsigned int> &neighbors) {
    neighbors.clear();
    std::vector<float> queryDistances{};
    queryDistances.resize(datasetSize);

    // get the exact NN of Q
    unsigned int index1;
    float distance_Q1 = HUGE_VAL;
    for (unsigned int index = 0; index < datasetSize; index++) {
        float const distance = sparseMatrix._computeDistance(queryIndex, index);
        queryDistances[index] = distance;
        if (index == queryIndex) continue;
        if (distance < distance_Q1) {
            distance_Q1 = distance;
            index1 = index;
        }
    }

    // check for interference by first neighbor NN
    neighbors.push_back(index1);
    std::vector<unsigned int> L1{};
    L1.reserve(datasetSize);  // single allocation for speed
    for (unsigned int index2 = 0; index2 < datasetSize; index2++) {
        if (index2 == queryIndex || index2 == index1) continue;
        float const distance_Q2 = queryDistances[index2];
        float const distance_12 = sparseMatrix._computeDistance(index1, index2);

        // add if not HSP inequalities NOT satisfied
        if (distance_Q1 >= distance_Q2 || distance_12 >= distance_Q2) {
            L1.push_back(index2);
        }
    }

    // now, find all HSP neighbors remaining
    unsigned int M = (unsigned int)L1.size();
    while (M > 0) {
        // find next nearest neighbor
        unsigned int index1;
        float distance_Q1 = HUGE_VAL;
        for (unsigned int i1 = 0; i1 < M; i1++) {
            unsigned int index = L1[i1];
            float const distance = queryDistances[index];
            if (distance < distance_Q1) {
                distance_Q1 = distance;
                index1 = index;
            }
        }
        neighbors.push_back(index1);

        // remove all that are interfered with
        std::vector<unsigned int>::iterator it2;
        for (it2 = L1.begin(); it2 != L1.end(); /* iterate in loop*/) {
            unsigned int index2 = (*it2);
            if (index2 == index1) {
                it2 = L1.erase(it2);  // remove hsp neighbor from list
                continue;
            }
            float const distance_12 = sparseMatrix._computeDistance(index1, index2);
            float const distance_Q2 = queryDistances[index2];

            // delete if HSP inequalities are satisfied
            if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                it2 = L1.erase(it2);
            } else {
                ++it2;
            }
        }
        M = (unsigned int)L1.size();
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

// compute distance using sparse matrix
float const GHSP::computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix) {
    return sparseMatrix._computeDistance(index1, index2);
}

// get or retrieve query distnace
float const GHSP::getQueryDistance(unsigned int const queryIndex, unsigned int const index2, SparseMatrix &sparseMatrix,
                                   std::vector<float> &distanceList) {
    float distance = distanceList[index2];
    if (distance < 0) {
        distance = sparseMatrix._computeDistance(queryIndex, index2);
        distanceList[index2] = distance;
    }
    return distance;
}

void GHSP::printSet(std::vector<unsigned int> const &set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}