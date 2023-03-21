#include "GHSP.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

// compute distance using sparse matrix
float const computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);

// get or retrieve query distnace
float const getQueryDistance(unsigned int const queryIndex, unsigned int const index2, SparseMatrix &sparseMatrix,
                             std::vector<float> &distanceList);

/**
 * @brief Check if a pivot is invalidated/safe against all existing HSP neighbors
 *
 * @param pivot             // the pivot being validated
 * @param queryIndex
 * @param neighbors         // HSP neighbors of queryIndex
 * @param sparseMatrix
 * @param queryDistances
 * @param pointValid        // is the point (index of pivot) itself invalidated? for NN update
 * @return 0: fully safe (not-intermediate), 1: invalidated (remove), 2: intermediate
 */
bool validatePoint(unsigned int const pointIndex, unsigned int const queryIndex,
                   std::vector<unsigned int> const &neighbors, SparseMatrix &sparseMatrix,
                   std::vector<float> &queryDistances) {
    float const distance_Q2 = getQueryDistance(queryIndex, pointIndex, sparseMatrix, queryDistances);

    // pointValid:
    //  true: the index of pivot is not invalidated by any HSP neighbor
    //  false: otherwise
    bool pointValid = true;

    // iterate through each HSP neighbor, check against GHSP inequalities
    for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
        unsigned int const x1_bar = neighbors[it1];
        float const distance_Q1 = getQueryDistance(queryIndex, x1_bar, sparseMatrix, queryDistances);
        float const distance_12 = computeDistance(x1_bar, pointIndex, sparseMatrix);

        // check point against normal HSP inequalities
        if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
            pointValid = false;
            break;
        }
    }

    return pointValid;
}

/**
 * @brief Check if a pivot is invalidated/safe against all existing HSP neighbors
 *
 * @param pivot             // the pivot being validated
 * @param queryIndex
 * @param neighbors         // HSP neighbors of queryIndex
 * @param sparseMatrix
 * @param queryDistances
 * @param pointValid        // is the point (index of pivot) itself invalidated? for NN update
 * @return 0: fully safe (not-intermediate), 1: invalidated (remove), 2: intermediate
 */
int validatePivot(Pivot const &pivot, unsigned int const queryIndex, std::vector<unsigned int> const &neighbors,
                  SparseMatrix &sparseMatrix, std::vector<float> &queryDistances, bool &pointValid) {
    float const distance_Q2 = getQueryDistance(queryIndex, pivot._index, sparseMatrix, queryDistances);
    float const radius = pivot._maxChildDistance;

    // pointValid:
    //  true: the index of pivot is not invalidated by any HSP neighbor
    //  false: otherwise
    pointValid = true;

    // flag_outcome:
    //      0: pivot domain fully safe by ALL neighbors, not intermediate
    //      1: pivot domain fully invalidated by a neighbor, remove from consideration
    //      2: neither, remains intermediate
    int flag_outcome = 0;

    // iterate through each HSP neighbor, check against GHSP inequalities
    for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
        unsigned int const x1_bar = neighbors[it1];
        float const distance_Q1 = getQueryDistance(queryIndex, x1_bar, sparseMatrix, queryDistances);
        float const distance_12 = computeDistance(x1_bar, pivot._index, sparseMatrix);

        // check point against normal HSP inequalities
        if (pointValid) {
            if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                pointValid = false;
            }
        }

        // check if the HSP neighbor fully invalidates the domain
        if (distance_Q1 < distance_Q2 - radius) {
            if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                flag_outcome = 1;
                break;
            }
        }

        // if still has chance to be fully safe
        if (flag_outcome == 0) {
            if (distance_Q1 > distance_Q2 + radius ||
                distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                continue;
            } else {
                flag_outcome = 2;
            }
        }
    }
    return flag_outcome;
}

/**
 * @brief Search to find the closest point to Q, recursive.
 * Only for non-intermediate points, does not employ validation checks
 *
 * @param parentPivot
 * @param queryIndex
 * @param index1
 * @param dmin
 * @param sparseMatrix
 * @param queryDistances
 */
void recursiveSearch(Pivot const &parentPivot, unsigned int const queryIndex,
                     std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
                     SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {
    if (parentPivot._childCount <= 0) return;

    // iterate through pivot domain
    std::vector<Pivot>::const_iterator it1;
    for (it1 = parentPivot._pivotDomain.begin(); it1 != parentPivot._pivotDomain.end(); it1++) {
        Pivot const &childPivot = (*it1);
        float const distance_Qc = getQueryDistance(queryIndex, childPivot._index, sparseMatrix, queryDistances);

        // update active NN only if its not an HSP neighbor already
        if (distance_Qc < dmin) {
            if (std::find(neighbors.begin(), neighbors.end(), childPivot._index) == neighbors.end()) {
                dmin = distance_Qc;
                index1 = childPivot._index;
            }
        }

        // check if domain can hold a closer point
        if (distance_Qc <= dmin + childPivot._maxChildDistance) {
            if (childPivot._childCount > 0) {
                recursiveSearch(childPivot, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances);
            }
        }
    }

    return;
}

/**
 * @brief Search to find the closest point to Q, recursive.
 * For intermediate pivots: requires a series of validation checks
 *
 * @param parentPivot
 * @param queryIndex
 * @param index1
 * @param dmin
 * @param sparseMatrix
 * @param queryDistances
 */
void recursiveSearch_Validation(Pivot const &parentPivot, unsigned int const queryIndex,
                                std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
                                SparseMatrix &sparseMatrix, std::vector<float> &queryDistances,
                                std::vector<Pivot> &newPivotList) {
    // iterate through pivot domain
    if (parentPivot._childCount <= 0) return;

    std::vector<Pivot>::const_iterator it1;
    for (it1 = parentPivot._pivotDomain.begin(); it1 != parentPivot._pivotDomain.end(); it1++) {
        Pivot childPivot = (*it1);
        float const distance_Qc = getQueryDistance(queryIndex, childPivot._index, sparseMatrix, queryDistances);

        // ensure if pivot/pivot domain is active/invalidated
        bool pointValid;
        int flag_outcome = validatePivot(childPivot, queryIndex, neighbors, sparseMatrix, queryDistances, pointValid);
        if (flag_outcome == 0) {  // entire domain safe
            childPivot.setIntermediate(false);
        } else if (flag_outcome == 1) {  // entire domain pivot is eliminated!
            continue;
        } else if (flag_outcome == 2) {  // pivot in intermediate region
            childPivot.setIntermediate(true);
        }

        // update dmin if active
        if (pointValid) {
            if (distance_Qc < dmin) {
                // update only if its not an HSP neighbor already
                if (std::find(neighbors.begin(), neighbors.end(), childPivot._index) == neighbors.end()) {
                    dmin = distance_Qc;
                    index1 = childPivot._index;
                }
            }
        }

        // check if domain can hold a closer point
        if (distance_Qc <= dmin + childPivot._maxChildDistance) {
            if (childPivot._childCount > 0) {
                if (flag_outcome == 0) {  // if active, don't need validation checks
                    newPivotList.push_back(childPivot);
                    recursiveSearch(childPivot, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances);

                } else if (flag_outcome == 2) {  // if intermeidate, need validation checks
                    recursiveSearch_Validation(childPivot, queryIndex, neighbors, index1, dmin, sparseMatrix,
                                               queryDistances, newPivotList);
                }
            } else {
                newPivotList.push_back(childPivot);
            }
        } else {
            // if we don't have to search the domain, add it to end of list (if int. or not)
            newPivotList.push_back(childPivot);
        }
    }

    return;
}

/**
 * @brief Performing GHSP search on a recursively-defined pivot index
 *
 * @param queryIndex
 * @param pivotsList
 * @param sparseMatrix
 * @param neighbors
 */
void GHSP::Recursive_GHSP_Search(unsigned int const queryIndex, std::vector<Pivot> const &pivotsList,
                                 SparseMatrix &sparseMatrix, std::vector<unsigned int> &neighbors) {
    neighbors.clear();

    // initialize query distance storage
    unsigned int datasetSize = sparseMatrix._datasetSize;
    std::vector<float> queryDistances{};
    queryDistances.resize(datasetSize, -1.0f);

    // one list to rule them all!
    std::vector<Pivot> A = pivotsList;

    // find the HSP neighbors in a loop
    while (A.size() > 0) {
        // TEMP: TIMING
        std::chrono::high_resolution_clock::time_point tStart, tEnd;
        tStart = std::chrono::high_resolution_clock::now();
        /**
         * ============================================================
         * @brief STEP 1: Collect dmin as preliminary closest active point
         * ============================================================
         */
        unsigned int index1;  // next hsp neighbor
        float dmin = HUGE_VAL;

        // update dmin with only validated pivots in I1
        std::vector<Pivot>::iterator it1;
        for (it1 = A.begin(); it1 != A.end(); it1++) {
            Pivot const &p = (*it1);
            float const distance_Qp = getQueryDistance(queryIndex, p._index, sparseMatrix, queryDistances);

            // if p is an hsp neighbor, ignore
            if (std::find(neighbors.begin(), neighbors.end(), p._index) != neighbors.end()) continue;

            // can update dmin if validated
            if (distance_Qp < dmin) {
                // validate against existing hsp neighbors
                bool flag_valid = true;
                if (p._intermediate) {
                    flag_valid = validatePoint(p._index, queryIndex, neighbors, sparseMatrix, queryDistances);
                }

                // update dmin if p2 is valid
                if (flag_valid) {
                    dmin = distance_Qp;
                    index1 = p._index;
                }
            }
        }

        /**
         * ============================================================
         * @brief STEP 2: Search the domains of pivots that could contain nearest neighbor
         * ============================================================
         */

        // Find the nearest neighbor by depth-first search of domains
        // Intermediate pivots will have to have their entire domain validated, thus such pivots are removed from the
        // list and their domain is added. We remain depth-first, so let us keep a list of new pivots we will insert at
        // the end.
        std::vector<Pivot> newPivotList{};
        for (it1 = A.begin(); it1 != A.end(); /* iterate in loop */) {
            Pivot const &p = (*it1);
            float const distance_Qp = getQueryDistance(queryIndex, p._index, sparseMatrix, queryDistances);

            // could pivot domain contain a closer point?
            if (distance_Qp > dmin + p._maxChildDistance) {  // No!
                ++it1;

            } else {  // yes, yes it can

                // need to recursively check all domains
                if (!p._intermediate) {  // no validation checks necessary
                    recursiveSearch(p, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances);
                    ++it1;

                } else {  // need validation checks!!
                    recursiveSearch_Validation(p, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances,
                                               newPivotList);

                    // since we had to validate the entire domain, remove this pivot from list
                    it1 = A.erase(it1);
                }
            }
        }

        // add all new pivots to the end
        A.insert(A.end(), newPivotList.begin(), newPivotList.end());
        if (dmin > 10000) continue;  // large, magic number signaling no neighbor possible

        // TEMP: TIMING
        tEnd = std::chrono::high_resolution_clock::now();
        double time_nns = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();
        printf("        - Time Search (ms): %.4f \n", time_nns * 1000);
        tStart = std::chrono::high_resolution_clock::now();

        // add the next HSP neighbor
        unsigned int const x1 = index1;
        float const distance_Q1 = dmin;
        neighbors.push_back(x1);

        printf("        * Length A: %u\n", A.size());

        /**
         * ============================================================
         * @brief STEP 3: Perform GHSP test on Active Pivots List
         * ============================================================
         */
        for (it1 = A.begin(); it1 != A.end(); /* iterate in loop */) {
            Pivot &p2 = (*it1);
            float const radius = p2._maxChildDistance;
            float const distance_Q2 = getQueryDistance(queryIndex, p2._index, sparseMatrix, queryDistances);
            float const distance_12 = computeDistance(x1, p2._index, sparseMatrix);

            // WHY DOES THIS BREAK IT
            if (radius <= 0) {  // normal HSP test
                if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                    it1 = A.erase(it1);
                } else {
                    ++it1;
                }
            } else {
                // Case 1: Entire domain invalidated
                if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1 &&
                    distance_Q1 < distance_Q2 - radius) {
                    it1 = A.erase(it1);
                }
                // Case 2: Entire domain safe
                else if (distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                    ++it1;
                }
                // Case 3: Intermediate pivot
                else {
                    p2.setIntermediate(true);
                    ++it1;
                }
            }
        }

        // TEMP: TIMING
        tEnd = std::chrono::high_resolution_clock::now();
        double time_val = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();
        printf("        - Time Val (ms): %.4f \n", time_val * 1000);
    }

    return;
}

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
        float const distance_Qp = getQueryDistance(queryIndex, pivot._index, sparseMatrix, queryDistances);

        // condition where pivot may contain nn
        if (distance_Qp <= searchRadius + radius) {
            std::vector<Pivot> const &pivotDomain = pivot._pivotDomain;
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

void GHSP::printSet(std::vector<unsigned int> const &set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}

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

// compute distance using sparse matrix
float const computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix) {
    return sparseMatrix._computeDistance(index1, index2);
}

// get or retrieve query distnace
float const getQueryDistance(unsigned int const queryIndex, unsigned int const index2, SparseMatrix &sparseMatrix,
                             std::vector<float> &distanceList) {
    float distance = distanceList[index2];
    if (distance < 0) {
        distance = sparseMatrix._computeDistance(queryIndex, index2);
        distanceList[index2] = distance;
    }
    return distance;
}