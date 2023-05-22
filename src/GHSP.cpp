#include "GHSP.hpp"

#include <cmath>
#include <queue>

namespace GHSP {
//> validation
bool validatePoint(unsigned int const pointIndex, unsigned int const queryIndex,
                   std::vector<unsigned int> const &neighbors, SparseMatrix &sparseMatrix);
int validatePivot(unsigned int const pivotIndex, float const radius, unsigned int const queryIndex,
                  std::vector<unsigned int> const &neighbors, SparseMatrix &sparseMatrix, bool &pointValid);

//> priority queue
typedef std::priority_queue<std::pair<float, unsigned int>, std::vector<std::pair<float, unsigned int>>,
                            std::greater<std::pair<float, unsigned int>>>
    PriorityQueue;

//> helper functions
float const computeDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);
float const getDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix);
;
void printSet(std::vector<unsigned int> const &set);
};  // namespace GHSP

/**
 * @brief 3-Layer GHSP Search
 *
 * @param queryIndex
 * @param pivotLayer
 * @param sparseMatrix
 * @param neighbors
 */
void GHSP::GHSP_3L(unsigned int const queryIndex, std::vector<PivotLayer> &pivotLayers, SparseMatrix &sparseMatrix,
                   std::vector<unsigned int> &neighbors) {
    neighbors.clear();

    // initialize query distance storage
    unsigned int datasetSize = sparseMatrix._datasetSize;
    sparseMatrix._clear();                      // remove all stored distances
    sparseMatrix._addNewReference(queryIndex);  // add vector for query distance storage

    // everything that can potentially become an HSP neighbor
    tsl::sparse_set<unsigned int> const &pivotIndices = *(pivotLayers[0].get_pivotIndices_ptr());
    PriorityQueue A1, I1, A2, I2, A3, I3{};
    PriorityQueue A1_copy, I1_copy, A2_copy, I2_copy, A3_copy, I3_copy{};
    PriorityQueue t2;
    std::pair<float, unsigned int> topOfQueue;
    int layerIndex = 0;

    // compute distance to all top layer pivots to initialize the queue
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        float const distance = getDistance(queryIndex, pivotIndex, sparseMatrix);
        A1.push(std::make_pair(distance, pivotIndex));
    }

    // Begin the GHSP Algorithm Loop
    int count = 0;
    while (A1.size() > 0 || I1.size() > 0 || A2.size() > 0 || I2.size() > 0 || A3.size() > 0) {
        /**
         * ============================================================
         * @brief STEP 1: Collect dmin as preliminary closest active point
         * ============================================================
         */

        //> Next HSP neighbor
        unsigned int index1;
        float dmin = HUGE_VAL;

        //> Update Dmin: Closest Active Point in L=3
        //------------------------------------------
        A3_copy = A3;
        while (A3_copy.size() > 0) {
            topOfQueue = A3_copy.top();
            A3_copy.pop();
            unsigned int const x2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;
                dmin = distance_Q2;
                index1 = x2;
            } else {
                break;
            }
        }

        //> Update Dmin: Closest Active Pivot in L=2
        //------------------------------------------
        A2_copy = A2;
        while (A2_copy.size() > 0) {
            topOfQueue = A2_copy.top();
            A2_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), p2) != neighbors.end()) continue;
                dmin = distance_Q2;
                index1 = p2;
            } else {
                break;
            }
        }

        //> Update Dmin: Closest Active Pivot in L=1
        //------------------------------------------
        A1_copy = A1;
        while (A1_copy.size() > 0) {
            topOfQueue = A1_copy.top();
            A1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), p2) != neighbors.end()) continue;
                dmin = distance_Q2;
                index1 = p2;
            } else {
                break;
            }
        }

        //> Update Dmin: Closest Intermediate Pivot in L=2
        //------------------------------------------------
        I2_copy = I2;
        while (I2_copy.size() > 0) {
            topOfQueue = I2_copy.top();
            I2_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            // can update dmin if validated
            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), p2) != neighbors.end()) continue;

                // update dmin if p2 is valid
                bool flag_validated = validatePoint(p2, queryIndex, neighbors, sparseMatrix);
                if (flag_validated) {
                    dmin = distance_Q2;
                    index1 = p2;
                }
            } else {  // no point in list closer
                break;
            }
        }

        //> Update Dmin: Closest Intermediate Pivot in L=1
        //------------------------------------------------
        I1_copy = I1;
        while (I1_copy.size() > 0) {
            topOfQueue = I1_copy.top();
            I1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            // can update dmin if validated
            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), p2) != neighbors.end()) continue;

                // update dmin if p2 is valid
                bool flag_validated = validatePoint(p2, queryIndex, neighbors, sparseMatrix);
                if (flag_validated) {
                    dmin = distance_Q2;
                    index1 = p2;
                }
            } else {  // no point in list closer
                break;
            }
        }

        /**
         * ============================================================================
         * @brief STEP 2: Search the domains of pivots that could contain next neighbor
         * ============================================================================
         */

        //> Search Domains: Active Pivots in L=2
        //--------------------------------------
        A2_copy = A2;
        while (A2_copy.size() > 0) {
            topOfQueue = A2_copy.top();
            A2_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius2 = pivotLayers[1].get_maxChildDistance(p2);

            // can it contain a closer point?
            if (distance_Q2 <= dmin + radius2) {
                std::vector<unsigned int> const &pivotDomain = pivotLayers[1].get_pivotChildren(p2);

                // now iterating on bottom layer
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const x2 = (*it3);
                    float const distance_Qx2 = getDistance(queryIndex, x2, sparseMatrix);

                    // this can update dmin!
                    if (distance_Qx2 < dmin) {
                        if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;
                        dmin = distance_Qx2;
                        index1 = x2;
                    }
                }
            }
        }

        //> Search Domains: Active Pivots in L=1
        //----------------------------------------------
        //  - first: iterate through top layer pivots and collect relevant domains
        t2 = PriorityQueue();  // temp list of active pivots in layer 2
        A1_copy = A1;
        while (A1_copy.size() > 0) {
            topOfQueue = A1_copy.top();
            A1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Qp2 = topOfQueue.first;
            float const radius1 = pivotLayers[0].get_maxChildDistance(p2);

            // Can Layer 1 Pivot Contain Closer Points?
            if (distance_Qp2 <= dmin + radius1) {
                std::vector<unsigned int> const &pivotDomain = pivotLayers[0].get_pivotChildren(p2);

                // now iterating through second layer pivots
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const p2c = (*it3);
                    float const distance_Q2c = getDistance(queryIndex, p2c, sparseMatrix);
                    float const radius2 = pivotLayers[1].get_maxChildDistance(p2c);

                    // this can update dmin!
                    if (distance_Q2c < dmin) {
                        if (std::find(neighbors.begin(), neighbors.end(), p2c) != neighbors.end()) continue;
                        dmin = distance_Q2c;
                        index1 = p2c;
                    }

                    // Add if it may contain a closer point
                    if (distance_Q2c <= dmin + radius2) {
                        t2.push(std::make_pair(distance_Q2c, p2c));
                    }
                }
            }
        }
        // - second: iterate through second layer pivots, searching relevant domains
        while (t2.size() > 0) {
            topOfQueue = t2.top();
            t2.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Qp2 = topOfQueue.first;
            float const radius2 = pivotLayers[1].get_maxChildDistance(p2);

            // Can Layer 1 Pivot Contain Closer Points?
            if (distance_Qp2 <= dmin + radius2) {
                std::vector<unsigned int> const &pivotDomain = pivotLayers[1].get_pivotChildren(p2);

                // now iterating on the bottom layer
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const x2 = (*it3);
                    float const distance_Qx2 = getDistance(queryIndex, x2, sparseMatrix);

                    // this can update dmin!
                    if (distance_Qx2 < dmin) {
                        if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;
                        dmin = distance_Qx2;
                        index1 = x2;
                    }
                }
            }
        }

        //> Search Domains: Intermediate Pivots in L=1
        //--------------------------------------------
        //  - first: iterate through intermediate top layer pivots and collect relevant domains
        I1_copy = I1;
        I1 = PriorityQueue();  // rebuild the list
        t2 = PriorityQueue();  // temporary list of active L2 pivots to search
        while (I1_copy.size() > 0) {
            topOfQueue = I1_copy.top();
            I1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius1 = pivotLayers[0].get_maxChildDistance(p2);

            // could pivot domain contain a closer point?
            if (distance_Q2 > dmin + radius1) {  // no! add back to I1
                I1.push(std::make_pair(distance_Q2, p2));

            } else {  // yes! validate entire domain, don't add back to I1
                std::vector<unsigned int> const &pivotDomain = pivotLayers[0].get_pivotChildren(p2);

                // iterating on second layer now
                std::vector<unsigned int>::const_iterator it2;
                for (it2 = pivotDomain.begin(); it2 != pivotDomain.end(); it2++) {
                    unsigned int const p2c = (*it2);
                    float const distance_Q2c = getDistance(queryIndex, p2c, sparseMatrix);
                    float const radius2 = pivotLayers[1].get_maxChildDistance(p2c);

                    // since not invalidated, add to active points list
                    bool flag_validated;
                    int flag_outcome = validatePivot(p2c, radius2, queryIndex, neighbors, sparseMatrix, flag_validated);
                    if (flag_validated) {
                        if (distance_Q2c < dmin) {
                            if (std::find(neighbors.begin(), neighbors.end(), p2c) == neighbors.end()) {
                                dmin = distance_Q2c;
                                index1 = p2c;
                            }
                        }
                    }

                    // based on entire pivot domain validation
                    if (flag_outcome == 0) {  // safe pivot
                        A2.push(std::make_pair(distance_Q2c, p2c));
                        if (distance_Q2c <= dmin + radius2) {  // YES
                            t2.push(std::make_pair(distance_Q2c, p2c));
                        }

                    } else if (flag_outcome == 1) {  // pivot and domain eliminated
                        continue;

                    } else if (flag_outcome == 2) {  // intermediate pivot and domain
                        I2.push(std::make_pair(distance_Q2c, p2c));
                    }
                }
            }
        }
        //  - second: iterate through new active 2nd layer pivots
        while (t2.size() > 0) {
            topOfQueue = t2.top();
            t2.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Qp2 = topOfQueue.first;
            float const radius2 = pivotLayers[1].get_maxChildDistance(p2);

            // Can Layer 3 Pivot Contain Closer Points?
            if (distance_Qp2 <= dmin + radius2) {
                std::vector<unsigned int> const &pivotDomain = pivotLayers[1].get_pivotChildren(p2);

                // now iterating on the bottom layer
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const x2 = (*it3);
                    float const distance_Qx2 = getDistance(queryIndex, x2, sparseMatrix);

                    // this can update dmin!
                    if (distance_Qx2 < dmin) {
                        if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;
                        dmin = distance_Qx2;
                        index1 = x2;
                    }
                }
            }
        }
        // - third: iterate through intermediate second level pivots
        I2_copy = I2;
        I2 = PriorityQueue();
        while (I2_copy.size() > 0) {
            topOfQueue = I2_copy.top();
            I2_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Qp2 = topOfQueue.first;
            float const radius2 = pivotLayers[1].get_maxChildDistance(p2);

            // Can Layer 2 Pivot Contain Closer Points?
            if (distance_Qp2 > dmin + radius2) {  // no!
                I2.push(std::make_pair(distance_Qp2, p2));

            } else {  // YES
                std::vector<unsigned int> const &pivotDomain = pivotLayers[1].get_pivotChildren(p2);

                // now iterating on the bottom layer
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const x2 = (*it3);
                    float const distance_Qx2 = getDistance(queryIndex, x2, sparseMatrix);

                    // update dmin if p2 is valid
                    bool flag_validated = validatePoint(x2, queryIndex, neighbors, sparseMatrix);
                    if (flag_validated) {
                        A3.push(std::make_pair(distance_Qx2, x2));

                        if (distance_Qx2 < dmin) {
                            if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;
                            dmin = distance_Qx2;
                            index1 = x2;
                        }
                    }
                }
            }
        }

        if (dmin > 10000) break;  // large, magic number signaling no new neighbor possible
        unsigned int const x1 = index1;
        float const distance_Q1 = getDistance(queryIndex, x1, sparseMatrix);
        neighbors.push_back(x1);
        sparseMatrix._addNewReference(x1);

        /**
         * ============================================================
         * @brief STEP 3: Perform GHSP test on all pivots agasinst GHSP between Q,x1
         * ============================================================
         */

        //> GHSP Test: Against Top Layer Intermediate Pivots
        //--------------------------------------------------
        // Three Options:
        //  1. Entire Pivot Domain Invalidated
        //  2. Pivot Remains in List
        //  3. Entire Domain must be Examined by Brute Force
        I1_copy = I1;
        I1 = PriorityQueue();
        while (I1_copy.size() > 0) {
            topOfQueue = I1_copy.top();
            I1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius = pivotLayers[0].get_maxChildDistance(p2);

            //> Pivot Domain cannot be excluded by x1 or any future HSP neighbor
            if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                I1.push(std::make_pair(distance_Q2, p2));  // TEMP
            }

            //> Otherwise, Check for Invalidation
            else {
                float const distance_12 = getDistance(x1, p2, sparseMatrix);

                // Proposition: Entire Domain Invalidated
                if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                    continue;
                }
                // Else, Remains in List
                else {
                    I1.push(std::make_pair(distance_Q2, p2));
                }
            }
        }

        //> GHSP Test: Against Second Layer Intermediate Pivots
        //-----------------------------------------------------
        // Three Options:
        //  1. Entire Pivot Domain Invalidated
        //  2. Pivot Remains in List
        //  3. Entire Domain must be Examined by Brute Force
        I2_copy = I2;
        I2 = PriorityQueue();
        while (I2_copy.size() > 0) {
            topOfQueue = I2_copy.top();
            I2_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius = pivotLayers[1].get_maxChildDistance(p2);

            //> Pivot Domain cannot be excluded by x1 or any future HSP neighbor
            if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                I2.push(std::make_pair(distance_Q2, p2));  // TEMP
            }

            //> Otherwise, Check for Invalidation
            else {
                float const distance_12 = getDistance(x1, p2, sparseMatrix);

                // Proposition: Entire Domain Invalidated
                if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                    continue;
                }
                // Else, Remains in List
                else {
                    I2.push(std::make_pair(distance_Q2, p2));
                }
            }
        }

        //> GHSP Test: Against Top Layer Active Pivots
        //--------------------------------------------
        // Iterate through active pivots
        // Four Options:
        //  1. Entire Pivot Domain Invalidated
        //  2. Pivot Remains in List
        //  3. Pivot Moves to Intermediate Pivot List
        //  4. Entire Domain must be Examined by Brute Force
        A1_copy = A1;
        A1 = PriorityQueue();
        while (A1_copy.size() > 0) {
            topOfQueue = A1_copy.top();
            A1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius = pivotLayers[0].get_maxChildDistance(p2);

            // Pivot cannot be excluded by x1 or any future HSP neighbor
            if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                I1.push(std::make_pair(distance_Q2, p2));  // TEMP
            }

            //> Otherwise, Check for Invalidation
            else {
                float const distance_12 = getDistance(x1, p2, sparseMatrix);

                // Proposition: Entire Domain Invalidated
                if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                    continue;
                }
                // Proposition: Entire Domain Safe
                else if (distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                    A1.push(std::make_pair(distance_Q2, p2));
                }
                // Otherwise, Pivot Falls in Intermediate Region
                else {
                    I1.push(std::make_pair(distance_Q2, p2));
                }
            }
        }

        //> GHSP Test: Against Second Layer Active Pivots
        //-----------------------------------------------------
        // Iterate through active pivots
        // Four Options:
        //  1. Entire Pivot Domain Invalidated
        //  2. Pivot Remains in List
        //  3. Pivot Moves to Intermediate Pivot List
        //  4. Entire Domain must be Examined by Brute Force
        A2_copy = A2;
        A2 = PriorityQueue();
        while (A2_copy.size() > 0) {
            topOfQueue = A2_copy.top();
            A2_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius = pivotLayers[1].get_maxChildDistance(p2);

            // Pivot cannot be excluded by x1 or any future HSP neighbor
            if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                I2.push(std::make_pair(distance_Q2, p2));  // TEMP
            }

            //> Otherwise, Check for Invalidation
            else {
                float const distance_12 = getDistance(x1, p2, sparseMatrix);

                // Proposition: Entire Domain Invalidated
                if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                    continue;
                }
                // Proposition: Entire Domain Safe
                else if (distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                    A2.push(std::make_pair(distance_Q2, p2));
                }
                // Otherwise, Pivot Falls in Intermediate Region
                else {
                    I2.push(std::make_pair(distance_Q2, p2));
                }
            }
        }

        //> GHSP Test: Against Bottom Layer Active Points
        //-----------------------------------------------------
        // Iterate through active pivots
        // Two Options:
        //  1. Entire Pivot Domain Invalidated
        //  2. Pivot Remains in List
        A3_copy = A3;
        A3 = PriorityQueue();
        while (A3_copy.size() > 0) {
            topOfQueue = A3_copy.top();
            A3_copy.pop();
            unsigned int const x2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;

            // HSP test
            float const distance_12 = getDistance(x1, x2, sparseMatrix);
            if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                continue;
            } else {
                A3.push(std::make_pair(distance_Q2, x2));
            }
        }
    }

    return;
}

/**
 * @brief
 *
 * @param pivotIndex        // pivot to be tested
 * @param radius            // radius of the pivot
 * @param queryIndex        //
 * @param neighbors         // existing HSP neighbors of Q
 * @param sparseMatrix
 * @param pointValid        // flag for if point itself is invalidated
 * @return int |--> 0: safe, 1: invalidated, 2: intermediate
 */
int GHSP::validatePivot(unsigned int const pivotIndex, float const radius, unsigned int const queryIndex,
                        std::vector<unsigned int> const &neighbors, SparseMatrix &sparseMatrix, bool &pointValid) {
    float const distance_Q2 = getDistance(queryIndex, pivotIndex, sparseMatrix);

    // pointValid:
    //  true: the index of pivot is not invalidated by any HSP neighbor
    //  false: else
    pointValid = true;

    // flag_outcome:
    //      0: pivot domain fully safe by ALL neighbors, not intermediate
    //      1: pivot domain fully invalidated by a neighbor, remove from consideration
    //      2: neither, remains intermediate
    int flag_outcome = 0;

    // iterate through each HSP neighbor, check against GHSP inequalities
    for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
        unsigned int const x1_bar = neighbors[it1];
        float const distance_Q1 = getDistance(queryIndex, x1_bar, sparseMatrix);
        float const distance_12 = getDistance(x1_bar, pivotIndex, sparseMatrix);

        // check point against normal HSP inequalities
        if (pointValid) {  // if not yet invalidated
            if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                pointValid = false;
            }
        }

        // check if the HSP neighbor fully invalidates the domain
        if (distance_Q1 < distance_Q2 - radius) {
            if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                pointValid = false;
                flag_outcome = 1;
                break;
            }
        }

        // check if the domain is safe
        if (flag_outcome == 0) {  // if not yet intermediate
            if (distance_Q1 > distance_Q2 + radius ||
                distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                continue;  // safe
            } else {
                flag_outcome = 2;
            }
        }
    }

    return flag_outcome;
}

/**
 * @brief 2-Layer GHSP Search
 *
 * @param queryIndex
 * @param pivotLayer
 * @param sparseMatrix
 * @param neighbors
 */
void GHSP::GHSP_2L(unsigned int const queryIndex, std::vector<PivotLayer> &pivotLayers, SparseMatrix &sparseMatrix,
                   std::vector<unsigned int> &neighbors) {
    neighbors.clear();
    PivotLayer& pivotLayer = pivotLayers[0];

    // initialize query distance storage
    unsigned int datasetSize = sparseMatrix._datasetSize;
    sparseMatrix._clear();                      // remove all stored distances
    sparseMatrix._addNewReference(queryIndex);  // add vector for query distance storage

    // everything that can potentially become an HSP neighbor
    tsl::sparse_set<unsigned int> const &pivotIndices = *pivotLayer.get_pivotIndices_ptr();
    PriorityQueue A1, I1, A2, I2{};
    PriorityQueue A1_copy, I1_copy, A2_copy, I2_copy{};
    std::pair<float, unsigned int> topOfQueue;

    // compute distance to all top layer pivots to initialize the queue
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        float const distance = getDistance(queryIndex, pivotIndex, sparseMatrix);
        A1.push(std::make_pair(distance, pivotIndex));
    }

    // Begin the GHSP Algorithm Loop
    while (A1.size() > 0 || I1.size() > 0 || A2.size() > 0) {
        /**
         * ============================================================
         * @brief STEP 1: Collect dmin as preliminary closest active point
         * ============================================================
         */

        //> Next HSP neighbor
        unsigned int index1;
        float dmin = HUGE_VAL;

        //> Update dmin with closest active point
        A2_copy = A2;
        while (A2_copy.size() > 0) {
            topOfQueue = A2_copy.top();
            A2_copy.pop();
            unsigned int const x2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;
                dmin = distance_Q2;
                index1 = x2;
            } else {
                break;
            }
        }

        //> Update dmin with closest active pivot
        A1_copy = A1;
        while (A1_copy.size() > 0) {
            topOfQueue = A1_copy.top();
            A1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), p2) != neighbors.end()) continue;
                dmin = distance_Q2;
                index1 = p2;
            } else {
                break;
            }
        }

        //> Update dmin with closest intermediate pivot
        I1_copy = I1;
        while (I1_copy.size() > 0) {
            topOfQueue = I1_copy.top();
            I1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;

            // can update dmin if validated
            if (distance_Q2 < dmin) {
                if (std::find(neighbors.begin(), neighbors.end(), p2) != neighbors.end()) continue;

                // Update dmin if p2 can be a potential neighbor (not invalidated)
                bool flag_validated = validatePoint(p2, queryIndex, neighbors, sparseMatrix);
                if (flag_validated) {
                    dmin = distance_Q2;
                    index1 = p2;
                }
            } else {  // no point in list closer
                break;
            }
        }

        /**
         * ============================================================
         * @brief STEP 2: Search the domains of pivots that could contain nearest neighbor
         * ============================================================
         */

        //> Update dmin with domains of active pivots
        A1_copy = A1;
        while (A1_copy.size() > 0) {
            topOfQueue = A1_copy.top();
            A1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Qp2 = topOfQueue.first;
            float const radius = pivotLayer.get_maxChildDistance(p2);

            if (distance_Qp2 <= dmin + radius) {
                std::vector<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);

                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const x2 = (*it3);
                    float const distance_Qx2 = getDistance(queryIndex, x2, sparseMatrix);

                    // this can update dmin!
                    if (distance_Qx2 < dmin) {
                        if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;
                        dmin = distance_Qx2;
                        index1 = x2;
                    }
                }
            }
        }

        // update dmin with domains of intermediate pivots
        I1_copy = I1;
        I1 = PriorityQueue();
        while (I1_copy.size() > 0) {
            topOfQueue = I1_copy.top();
            I1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Qp2 = topOfQueue.first;
            float const radius = pivotLayer.get_maxChildDistance(p2);

            // could pivot domain contain a closer point?
            if (distance_Qp2 > dmin + radius) {  // no! add back to I1
                I1.push(std::make_pair(distance_Qp2, p2));

            } else {  // yes! validate entire domain, don't add back to I1

                std::vector<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const x2 = (*it3);
                    float const distance_Q2 = getDistance(queryIndex, x2, sparseMatrix);
                    if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;

                    // validate point against all previous HSP neighbors
                    bool flag_validated = true;
                    for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
                        unsigned int const x1_bar = neighbors[it1];
                        float const distance_Q1 = getDistance(queryIndex, x1_bar, sparseMatrix);
                        float const distance_12 = getDistance(x1_bar, x2, sparseMatrix);
                        if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                            flag_validated = false;
                            break;
                        }
                    }

                    // since not invalidated, add to active points list
                    if (flag_validated) {
                        A2.push(std::make_pair(distance_Q2, x2));

                        // this can update dmin!
                        if (distance_Q2 < dmin) {
                            dmin = distance_Q2;
                            index1 = x2;
                        }
                    }
                }
            }
        }
        if (dmin > 10000) continue;  // large, magic number signaling no neighbor possible
        unsigned int const x1 = index1;
        float const distance_Q1 = getDistance(queryIndex, x1, sparseMatrix);
        neighbors.push_back(x1);
        sparseMatrix._addNewReference(x1);

        /**
         * ============================================================
         * @brief STEP 3: Perform GHSP test on all pivots agasinst GHSP between Q,x1
         * ============================================================
         */

        // Iterate through intermediate pivots
        // Three Options:
        //  1. Entire Pivot Domain Invalidated
        //  2. Pivot Remains in List
        //  3. Entire Domain must be Examined by Brute Force
        I1_copy = I1;
        I1 = PriorityQueue();
        I2 = PriorityQueue();
        while (I1_copy.size() > 0) {
            topOfQueue = I1_copy.top();
            I1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius = pivotLayer.get_maxChildDistance(p2);

            //> Pivot Domain cannot be excluded by x1 or any future HSP neighbor
            if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                std::vector<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);

                // Entire domain must be validated against all prior HSP neighbor
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const pointIndex = (*it3);
                    float const distance = getDistance(queryIndex, pointIndex, sparseMatrix);
                    I2.push(std::make_pair(distance, pointIndex));
                }
            }

            //> Otherwise, Check for Invalidation
            else {
                float const distance_12 = getDistance(x1, p2, sparseMatrix);

                // Proposition: Entire Domain Invalidated
                if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                    continue;
                }
                // Else, Remains in List
                else {
                    I1.push(std::make_pair(distance_Q2, p2));
                }
            }
        }

        // Iterate through active pivots
        // Four Options:
        //  1. Entire Pivot Domain Invalidated
        //  2. Pivot Remains in List
        //  3. Pivot Moves to Intermediate Pivot List
        //  4. Entire Domain must be Examined by Brute Force

        A1_copy = A1;
        A1 = PriorityQueue();
        while (A1_copy.size() > 0) {
            topOfQueue = A1_copy.top();
            A1_copy.pop();
            unsigned int const p2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            float const radius = pivotLayer.get_maxChildDistance(p2);

            // Pivot cannot be excluded by x1 or any future HSP neighbor
            if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                std::vector<unsigned int> const &pivotDomain = pivotLayer.get_pivotChildren(p2);

                // Entire domain must be validated against all prior HSP neighbor
                std::vector<unsigned int>::const_iterator it3;
                for (it3 = pivotDomain.begin(); it3 != pivotDomain.end(); it3++) {
                    unsigned int const pointIndex = (*it3);
                    float const distance = getDistance(queryIndex, pointIndex, sparseMatrix);
                    A2.push(std::make_pair(distance, pointIndex));
                }
            }

            //> Otherwise, Check for Invalidation
            else {
                float const distance_12 = getDistance(x1, p2, sparseMatrix);

                // Proposition: Entire Domain Invalidated
                if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                    continue;
                }
                // Proposition: Entire Domain Safe
                else if (distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                    A1.push(std::make_pair(distance_Q2, p2));
                }
                // Otherwise, Pivot Falls in Intermediate Region
                else {
                    I1.push(std::make_pair(distance_Q2, p2));
                }
            }
        }

        /**
         * ============================================================
         * @brief STEP 4: Perform normal HSP test on all points in with HSP(Q,x1)
         * ============================================================
         */

        A2_copy = A2;
        A2 = PriorityQueue();
        while (A2_copy.size() > 0) {
            topOfQueue = A2_copy.top();
            A2_copy.pop();
            unsigned int const x2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;

            // HSP test
            float const distance_12 = getDistance(x1, x2, sparseMatrix);
            if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                continue;
            } else {
                A2.push(std::make_pair(distance_Q2, x2));
            }
        }

        /**
         * ============================================================
         * @brief STEP 7: Validate all points in I2 against all HSP neighbors
         * ============================================================
         */

        I2_copy = I2;
        I2 = PriorityQueue();
        while (I2_copy.size() > 0) {
            topOfQueue = I2_copy.top();
            I2_copy.pop();
            unsigned int const x2 = topOfQueue.second;
            float const distance_Q2 = topOfQueue.first;
            if (std::find(neighbors.begin(), neighbors.end(), x2) != neighbors.end()) continue;

            // validate against all existing hsp neighbors
            bool flag_validated = true;
            for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
                unsigned int const x1_bar = neighbors[it1];
                float const distance_Q1 = getDistance(queryIndex, x1_bar, sparseMatrix);
                float const distance_12 = getDistance(x1_bar, x2, sparseMatrix);
                if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                    flag_validated = false;
                    break;
                }
            }

            if (flag_validated) {
                A2.push(std::make_pair(distance_Q2, x2));
            }
        }
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
                         std::vector<unsigned int> const &neighbors, SparseMatrix &sparseMatrix) {
    float const distance_Q2 = getDistance(queryIndex, pointIndex, sparseMatrix);

    // pointValid:
    //  true: the index of pivot is not invalidated by any HSP neighbor
    //  false: otherwise
    bool pointValid = true;

    // iterate through each HSP neighbor, check against GHSP inequalities
    for (unsigned int it1 = 0; it1 < neighbors.size(); it1++) {
        unsigned int const x1_bar = neighbors[it1];
        float const distance_Q1 = getDistance(queryIndex, x1_bar, sparseMatrix);
        float const distance_12 = getDistance(x1_bar, pointIndex, sparseMatrix);

        // check point against normal HSP inequalities
        if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
            pointValid = false;
            break;
        }
    }

    return pointValid;
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
void GHSP::HSP(unsigned int const queryIndex, unsigned int const datasetSize, SparseMatrix &sparseMatrix,
               std::vector<unsigned int> &neighbors) {
    neighbors.clear();
    sparseMatrix._clear();                      // remove all stored distances
    sparseMatrix._addNewReference(queryIndex);  // add vector for query distance storage

    // get the exact NN of Q
    unsigned int index1;
    float distance_Q1 = HUGE_VAL;
    for (unsigned int index = 0; index < datasetSize; index++) {
        float const distance = getDistance(queryIndex, index, sparseMatrix);
        if (index == queryIndex) continue;
        if (distance < distance_Q1) {
            distance_Q1 = distance;
            index1 = index;
        }
    }
    neighbors.push_back(index1);

    // check for interference by first neighbor NN
    std::vector<unsigned int> L1{};
    L1.reserve(datasetSize);  // single allocation for speed
    for (unsigned int index2 = 0; index2 < datasetSize; index2++) {
        if (index2 == queryIndex || index2 == index1) continue;
        float const distance_Q2 = getDistance(queryIndex, index2, sparseMatrix);
        float const distance_12 = computeDistance(index1, index2, sparseMatrix);

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
            float const distance = getDistance(queryIndex, index, sparseMatrix);
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
            float const distance_Q2 = getDistance(queryIndex, index2, sparseMatrix);
            float const distance_12 = computeDistance(index1, index2, sparseMatrix);

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
float const GHSP::getDistance(unsigned int const index1, unsigned int const index2, SparseMatrix &sparseMatrix) {
    return sparseMatrix._getDistance(index1, index2);
}

void GHSP::printSet(std::vector<unsigned int> const &set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}