#include "GHSP.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace GHSP {
bool validatePoint(unsigned int const pointIndex, unsigned int const queryIndex,
                   std::vector<unsigned int> const &neighbors, SparseMatrix &sparseMatrix,
                   std::vector<float> &queryDistances);
int validatePivot(Pivot const &pivot, unsigned int const queryIndex, std::vector<unsigned int> const &neighbors,
                  SparseMatrix &sparseMatrix, std::vector<float> &queryDistances, bool &pointValid);

void DepthFirstSearch(const Pivot* pivot, bool flag_intermediate, unsigned int const queryIndex, std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
                     SparseMatrix &sparseMatrix, std::vector<float> &queryDistances);

void BreadthFirstSearch(std::vector<const Pivot *> const &activePivots,
                           std::vector<const Pivot *> const &intermediatePivots,
                           std::vector<std::pair<const Pivot*, bool>> &pivotsToAdd, unsigned int const queryIndex,
                           std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
                           SparseMatrix &sparseMatrix, std::vector<float> &queryDistances);
}  // namespace GHSP

/**
 * @brief Performing GHSP search on a recursively-defined pivot index
 *
 * @param queryIndex
 * @param pivotsList
 * @param sparseMatrix
 * @param neighbors
 */
void GHSP::GHSP_Search(unsigned int const queryIndex, std::vector<Pivot> const &pivotsList, SparseMatrix &sparseMatrix,
                       std::vector<unsigned int> &neighbors) {
    neighbors.clear();

    // initialize query distance storage
    unsigned int datasetSize = sparseMatrix._datasetSize;
    std::vector<float> queryDistances{};
    queryDistances.resize(datasetSize, -1.0f);

    // one list to rule them all!
    // pair: pointer to pivot, flag for intermeidate or not.
    std::vector<std::pair<Pivot const *, bool>> A{};
    A.reserve(pivotsList.size());
    for (int it1 = 0; it1 < pivotsList.size(); it1++) {
        A.emplace_back(&pivotsList[it1], false);  // all start as not intermediate
    }

    // find the HSP neighbors in a loop
    std::vector<Pivot*> ActivePivotList{};
    while (A.size() > 0) {
        unsigned int index1;  // next hsp neighbor
        float dmin = HUGE_VAL;

        /**
         * =================================================================================
         * @brief STEP 1: Greedy Depth-First Search for Closest Active Point
         * =================================================================================
         */
        const Pivot* closestPivot;
        float closestPivotDistance = HUGE_VAL;
        bool closestPivotIntermediateFlag = true;
        std::vector<std::pair<const Pivot*, bool>>::iterator it1;
        for (it1 = A.begin(); it1 != A.end(); it1++) {
            Pivot const *pivot = (*it1).first;
            bool flag_intermediate = (*it1).second;
            float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

            if (distance < closestPivotDistance) {
                if (std::find(neighbors.begin(), neighbors.end(), pivot->_index) != neighbors.end()) continue;

                bool flag_valid = true;
                if (flag_intermediate) {  // validate against existing hsp neighbors
                    flag_valid = validatePoint(pivot->_index, queryIndex, neighbors, sparseMatrix, queryDistances);
                }

                // update dmin if p2 is valid
                if (flag_valid) {
                    closestPivotDistance = distance;
                    closestPivot = pivot;
                    closestPivotIntermediateFlag = flag_intermediate;

                    if (distance < dmin) {
                        dmin = distance;
                        index1 = pivot->_index;
                    }
                }
            }
        }

        // greedy depth-first search
        // ADJUST FOR INTERMEDIATE PIVOTS
        DepthFirstSearch(closestPivot, closestPivotIntermediateFlag, queryIndex, neighbors, index1, dmin,
                     sparseMatrix, queryDistances);

        /**
         * =================================================================================
         * @brief STEP 2: Breadth-First Search To Reduce Domains Examined
         * =================================================================================
         */

        // Find the nearest neighbor by breadth-first search of domains
        // Intermediate pivots will have to have their entire domain validated, thus such pivots are removed from the
        // list and their domain is added. 
        std::vector<const Pivot*> ActivePivotsList{};
        std::vector<const Pivot*> IntermediatePivotsList{};
        std::vector<std::pair<const Pivot*, bool>> newPivotList{};
        std::vector<Pivot>::const_iterator it2;
        for (it1 = A.begin(); it1 != A.end(); /* iterate in loop */) {
            const Pivot* pivot = (*it1).first;
            bool flag_intermediate = (*it1).second;
            float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

            // could pivot domain contain a closer point?
            if (distance <= dmin + pivot->_radius) {  // yes, yes it can

                // need to recursively check all domains
                if (!flag_intermediate) {  // no validation checks necessary
                    for (it2 = pivot->_pivotDomain.begin(); it2 != pivot->_pivotDomain.end(); it2++) {
                        ActivePivotsList.push_back(&(*it2));
                    }
                    ++it1;
                } else {  // need validation checks!!
                    for (it2 = pivot->_pivotDomain.begin(); it2 != pivot->_pivotDomain.end(); it2++) {
                        IntermediatePivotsList.push_back(&(*it2));
                    }
                    // since we have to validate the entire domain, remove this pivot from list
                    it1 = A.erase(it1);
                }
            } else {  // Yes
                ++it1;
            }
        }

        // Recursive Breadth-First Search
        BreadthFirstSearch(ActivePivotsList, IntermediatePivotsList, newPivotList, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances);
        A.insert(A.end(), newPivotList.begin(), newPivotList.end());
        if (dmin > 10000) continue;  // large, magic number signaling no neighbor possible

        /**
         * ========================================================================
         * @brief STEP 3: GHSP Test Against All Pivots with New HSP Neighbor
         * ========================================================================
         */
        unsigned int const x1 = index1;
        float const distance_Q1 = dmin;
        neighbors.push_back(x1);

        // iterate through each pivot in the list
        for (it1 = A.begin(); it1 != A.end(); /* iterate in loop */) {
            Pivot const &p2 = *((*it1).first);
            bool &flag_intermediate = (*it1).second;
            float const radius = p2._radius;
            float const distance_Q2 = getQueryDistance(queryIndex, p2._index, sparseMatrix, queryDistances);

            // call it intermediate and continue
            if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                flag_intermediate = true;
                ++it1;

            } else {
                float const distance_12 = computeDistance(x1, p2._index, sparseMatrix);

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
                    flag_intermediate = true;
                    ++it1;
                }
            }
        }
    }

    return;
}


/**
 * @brief Given pivots, breadth-first approach to NNS. Update dmin using all pivots. Collect select pivot domains for
 * search.
 *
 * @param activePivots                  // list of pivots in layer that are confirmed valid.
 * @param intermediatePivots            // list of pivots in layer requiring validation.
 * @param pivotsToAdd                   // list of pivots to add to ongoing list
 * @param queryIndex
 * @param neighbors
 * @param index1
 * @param dmin
 * @param sparseMatrix
 * @param queryDistances
 */
void GHSP::DepthFirstSearch(const Pivot* pivot, bool flag_intermediate, unsigned int const queryIndex, std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
                     SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {
    // store the closest point
    const Pivot* closestPivot;
    float closestPivotDistance = HUGE_VAL;

    std::vector<Pivot>::const_iterator it1;
    for (it1 = pivot->_pivotDomain.begin(); it1 != pivot->_pivotDomain.end(); it1++) {
        const Pivot* childPivot = &(*it1);
        float const distance = getQueryDistance(queryIndex, childPivot->_index, sparseMatrix, queryDistances);

        if (distance < closestPivotDistance) {
            if (std::find(neighbors.begin(), neighbors.end(), childPivot->_index) != neighbors.end()) continue;

            bool flag_valid = true;
            if (flag_intermediate) {  // validate against existing hsp neighbors
                flag_valid = validatePoint(childPivot->_index, queryIndex, neighbors, sparseMatrix, queryDistances);
            }

            // update dmin if p2 is valid
            if (flag_valid) {
                closestPivotDistance = distance;
                closestPivot = childPivot;

                // update closest pivot
                if (distance < dmin) {
                    dmin = distance;
                    index1 = childPivot->_index;
                }
            }
        }
    }

    // go down
    if (closestPivot->_childCount > 0) {
        DepthFirstSearch(closestPivot, flag_intermediate, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances);
    }

    return;
}




/**
 * @brief Given pivots, breadth-first approach to NNS. Update dmin using all pivots. Collect select pivot domains for
 * search.
 *
 * @param activePivots                  // list of pivots in layer that are confirmed valid.
 * @param intermediatePivots            // list of pivots in layer requiring validation.
 * @param pivotsToAdd                   // list of pivots to add to ongoing list
 * @param queryIndex
 * @param neighbors
 * @param index1
 * @param dmin
 * @param sparseMatrix
 * @param queryDistances
 */
void GHSP::BreadthFirstSearch(std::vector<const Pivot *> const &activePivots,
                           std::vector<const Pivot *> const &intermediatePivots,
                           std::vector<std::pair<const Pivot*, bool>> &pivotsToAdd, unsigned int const queryIndex,
                           std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
                           SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {
    // update dmin with list of active pivots
    std::vector<const Pivot *>::const_iterator it1;
    for (it1 = activePivots.begin(); it1 != activePivots.end(); it1++) {
        const Pivot *pivot = (*it1);
        float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

        // check for updating dmin
        if (std::find(neighbors.begin(), neighbors.end(), pivot->_index) != neighbors.end()) continue;
        if (distance < dmin) {
            dmin = distance;
            index1 = pivot->_index;
        }
    }

    // update dmin with list of intermediate pivots.
    std::vector<int> intermediatePivotStatus(intermediatePivots.size(),2);
    int it1_loc = 0;
    for (it1 = intermediatePivots.begin(); it1 != intermediatePivots.end(); it1++) {
        const Pivot *pivot = (*it1);
        float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

        // check status of pivot domain
        bool pointValid;
        bool flag_intermediate = true;
        int flag_outcome = validatePivot(*pivot, queryIndex, neighbors, sparseMatrix, queryDistances, pointValid);
        intermediatePivotStatus[it1_loc++] = flag_outcome;

        // can update dmin if validated
        if (pointValid) {
            if (distance < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), pivot->_index) == neighbors.end()) {
                    dmin = distance;
                    index1 = pivot->_index;
                }
            }
        }
    }

    // now, iterate through intermediate pivots: update dmin, collect domains
    it1_loc = 0;
    std::vector<Pivot>::const_iterator it2;
    std::vector<const Pivot *> newActivePivots{};
    std::vector<const Pivot *> newIntermediatePivots{};
    for (it1 = intermediatePivots.begin(); it1 != intermediatePivots.end(); it1++) {
        const Pivot *pivot = (*it1);
        float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

        // get the outcome from before
        int flag_outcome = intermediatePivotStatus[it1_loc++];
        bool flag_intermediate;
        if (flag_outcome == 0) {  // entire domain safe
            flag_intermediate = false;
        } else if (flag_outcome == 1) {  // entire domain pivot is eliminated!
            continue;
        } else if (flag_outcome == 2) {  // pivot in intermediate region
            flag_intermediate = true;
        }

        // could pivot domain contain a closer point?
        if (distance > dmin + pivot->_radius) {  // NO! All good
            pivotsToAdd.emplace_back(pivot, flag_intermediate);
            continue;
        };
        if (pivot->_childCount <= 0) {
            pivotsToAdd.emplace_back(pivot, flag_intermediate);
            continue;
        };

        // if not invalidated, add whole domain to active pivots. Add each to list.
        if (!flag_intermediate) {
            for (it2 = pivot->_pivotDomain.begin(); it2 != pivot->_pivotDomain.end(); it2++) {
                const Pivot *childPivot = &(*it2);
                newActivePivots.push_back(childPivot);
                pivotsToAdd.emplace_back(childPivot, false);
            }
        }
        // otherwise, add whole domain to intermeidate. domain must be examined, so don't add to list.
        else {
            for (it2 = pivot->_pivotDomain.begin(); it2 != pivot->_pivotDomain.end(); it2++) {
                newIntermediatePivots.push_back(&(*it2));
            }
        }
    }

    // finally, iterate through active pivots and collect domains that can contribute to dmin
    for (it1 = activePivots.begin(); it1 != activePivots.end(); it1++) {
        const Pivot *pivot = (*it1);
        float const distance = getQueryDistance(queryIndex, pivot->_index, sparseMatrix, queryDistances);

        // could pivot domain contain a closer point?
        if (distance <= dmin + pivot->_radius) {  // yes, yes it can
            if (pivot->_childCount > 0) {
                for (it2 = pivot->_pivotDomain.begin(); it2 != pivot->_pivotDomain.end(); it2++) {
                    newActivePivots.push_back(&(*it2));
                }
            }
        }
    }

    // Keep going down :)
    if (newActivePivots.size() > 0 || newIntermediatePivots.size() > 0) {
        BreadthFirstSearch(newActivePivots, newIntermediatePivots, pivotsToAdd, queryIndex, neighbors, index1, dmin,
                        sparseMatrix, queryDistances);
    }

    return;
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
bool GHSP::validatePoint(unsigned int const pointIndex, unsigned int const queryIndex,
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
int GHSP::validatePivot(Pivot const &pivot, unsigned int const queryIndex, std::vector<unsigned int> const &neighbors,
                        SparseMatrix &sparseMatrix, std::vector<float> &queryDistances, bool &pointValid) {
    float const distance_Q2 = getQueryDistance(queryIndex, pivot._index, sparseMatrix, queryDistances);
    float const radius = pivot._radius;  //_maxChildDistance;

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
 * ==============================================================================
 *
 *                      BRUTE FORCE
 *
 * ==============================================================================
 */

/**
 * @brief perform normal HSP Test
 *
 * @param queryIndex
 * @param datasetSize
 * @param sparseMatrix
 * @param neighbors
 */
void GHSP::HSP_Search(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
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










/**
 * @brief 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 *              OLD
 * 
 */

// /**
//  * @brief Search to find the closest point to Q, recursive.
//  * Only for non-intermediate points, does not employ validation checks
//  *
//  * @param parentPivot
//  * @param queryIndex
//  * @param index1
//  * @param dmin
//  * @param sparseMatrix
//  * @param queryDistances
//  */
// void GHSP::recursiveSearch(std::vector<const Pivot *> const &pivotList, unsigned int const queryIndex,
//                            std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
//                            SparseMatrix &sparseMatrix, std::vector<float> &queryDistances) {
//     // iterate through pivot domain
//     std::vector<Pivot>::const_iterator it1;
//     for (it1 = parentPivot._pivotDomain.begin(); it1 != parentPivot._pivotDomain.end(); it1++) {
//         Pivot const &childPivot = (*it1);
//         float const distance_Qc = getQueryDistance(queryIndex, childPivot._index, sparseMatrix, queryDistances);

//         // update active NN only if its not an HSP neighbor already
//         if (distance_Qc < dmin) {
//             if (std::find(neighbors.begin(), neighbors.end(), childPivot._index) == neighbors.end()) {
//                 dmin = distance_Qc;
//                 index1 = childPivot._index;
//             }
//         }

//         // check if domain can hold a closer point
//         if (distance_Qc <= dmin + childPivot._radius) {
//             if (childPivot._childCount > 0) {
//                 recursiveSearch(childPivot, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances);
//             }
//         }
//     }

//     return;
// }

// /**
//  * @brief Search to find the closest point to Q, recursive.
//  * For intermediate pivots: requires a series of validation checks
//  *
//  * @param parentPivot
//  * @param queryIndex
//  * @param index1
//  * @param dmin
//  * @param sparseMatrix
//  * @param queryDistances
//  */
// void GHSP::recursiveSearch_Validation(Pivot const &parentPivot, unsigned int const queryIndex,
//                                       std::vector<unsigned int> const &neighbors, unsigned int &index1, float &dmin,
//                                       SparseMatrix &sparseMatrix, std::vector<float> &queryDistances,
//                                       std::vector<std::pair<Pivot const *, bool>> &newPivotList) {
//     // iterate through pivot domain
//     // if (parentPivot._childCount <= 0) return;

//     std::vector<Pivot>::const_iterator it1;
//     for (it1 = parentPivot._pivotDomain.begin(); it1 != parentPivot._pivotDomain.end(); it1++) {
//         Pivot const &childPivot = (*it1);
//         bool flag_childIntermediate = false;
//         float const distance_Qc = getQueryDistance(queryIndex, childPivot._index, sparseMatrix, queryDistances);

//         // ensure if pivot/pivot domain is active/invalidated
//         bool pointValid;
//         int flag_outcome = validatePivot(childPivot, queryIndex, neighbors, sparseMatrix, queryDistances, pointValid);
//         if (flag_outcome == 0) {  // entire domain safe
//             flag_childIntermediate = false;
//         } else if (flag_outcome == 1) {  // entire domain pivot is eliminated!
//             continue;
//         } else if (flag_outcome == 2) {  // pivot in intermediate region
//             flag_childIntermediate = true;
//         }

//         // update dmin if active
//         if (pointValid) {
//             if (distance_Qc < dmin) {
//                 // update only if its not an HSP neighbor already
//                 if (std::find(neighbors.begin(), neighbors.end(), childPivot._index) == neighbors.end()) {
//                     dmin = distance_Qc;
//                     index1 = childPivot._index;
//                 }
//             }
//         }

//         // check if domain can hold a closer point
//         if (distance_Qc <= dmin + childPivot._radius) {
//             if (childPivot._childCount > 0) {
//                 if (!flag_childIntermediate) {  // if active, don't need validation checks
//                     newPivotList.emplace_back(&childPivot, flag_childIntermediate);
//                     recursiveSearch(childPivot, queryIndex, neighbors, index1, dmin, sparseMatrix, queryDistances);

//                 } else if (flag_childIntermediate) {  // if intermeidate, need validation checks
//                     recursiveSearch_Validation(childPivot, queryIndex, neighbors, index1, dmin, sparseMatrix,
//                                                queryDistances, newPivotList);
//                 }
//             } else {
//                 newPivotList.emplace_back(&childPivot, flag_childIntermediate);
//             }
//         } else {
//             // if we don't have to search the domain, add it to end of list (if int. or not)
//             newPivotList.emplace_back(&childPivot, flag_childIntermediate);
//         }
//     }

//     return;
// }