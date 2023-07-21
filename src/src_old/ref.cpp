void GHSP::CHSP(unsigned int const queryIndex, std::vector<PivotLayer> &pivotLayers, SparseMatrix &sparseMatrix,
                std::vector<unsigned int> &neighbors) {
    neighbors.clear();
    int const numLayers = pivotLayers.size() + 1;

    //  - initialize query distance storage
    unsigned int datasetSize = sparseMatrix._datasetSize;
    sparseMatrix._clear();                      // remove all stored distances
    sparseMatrix._addNewReference(queryIndex);  // add vector for query distance storage

    //  - initialize Lists for Multi-Layer HSP Search
    std::vector<PriorityQueue> activePivots{};  // pivots safe from invalidation by any existing neighbors
    activePivots.resize(numLayers);
    std::vector<PriorityQueue> intermediatePivots{};  // pivots questionable about safety
    intermediatePivots.resize(numLayers);

    //  - compute distance to all top layer pivots to initialize the queue
    tsl::sparse_set<unsigned int> const &pivotIndices = *(pivotLayers[0].get_pivotIndices_ptr());
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        float const distance = getDistance(queryIndex, pivotIndex, sparseMatrix);
        activePivots[0].push(std::make_pair(distance, pivotIndex));
    }

    //  - count the total number of points in lists
    unsigned int remainingPoints = 0;
    for (int layerIndex = 0; layerIndex < numLayers; layerIndex++) {
        remainingPoints += (activePivots[layerIndex].size() + intermediatePivots[layerIndex].size());
    }

    // printf("Made it...\n");

    // Begin the Cover Tree- HSP Algorithm Loop
    while (remainingPoints >= 0) {
        //
        /**
         * ======================================================
         *   STEP 1: Active Nearest Neighbor Search
         *     - find index1 as the next hsp neighbor of Q
         *
         * ======================================================
         */
        unsigned int index1;
        float dmin = HUGE_VAL;

        //> Check the active pivots bottom-up
        for (int layerIndex = numLayers - 1; layerIndex >= 0; layerIndex--) {
            // printf("- AP: Layer %d...\n",layerIndex);
            PriorityQueue list_copy = activePivots[layerIndex];

            // iterate through this copy of the queue
            std::pair<float, unsigned int> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                unsigned int const pivotIndex = topOfQueue.second;
                float const distance = topOfQueue.first;
                list_copy.pop();

                if (distance < dmin) {
                    if (std::find(neighbors.begin(), neighbors.end(), pivotIndex) != neighbors.end()) continue;
                    dmin = distance;
                    index1 = pivotIndex;
                } else {
                    break;
                }
            }
        }

        //> Check the intermediate pivots bottom-up
        for (int layerIndex = numLayers - 2; layerIndex >= 0; layerIndex--) {
            // printf("- IP: Layer %d...\n",layerIndex);
            PriorityQueue list_copy = intermediatePivots[layerIndex];

            // iterate through this copy of the queue
            std::pair<float, unsigned int> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                unsigned int const pivotIndex = topOfQueue.second;
                float const distance = topOfQueue.first;
                list_copy.pop();

                if (distance < dmin) {
                    if (std::find(neighbors.begin(), neighbors.end(), pivotIndex) != neighbors.end()) continue;

                    // update dmin if p2 is valid
                    bool flag_validated = validatePoint(pivotIndex, queryIndex, neighbors, sparseMatrix);
                    if (flag_validated) {
                        dmin = distance;
                        index1 = pivotIndex;
                    }
                } else {
                    break;
                }
            }
        }

        //> Check the Domains of Active Pivots Top-Down
        std::vector<PriorityQueue> copy_activePivots = activePivots;
        for (int layerIndex = 0; layerIndex < numLayers - 1; layerIndex++) {
            // printf("- APD: Layer %d...\n",layerIndex);
            //
            // iterate through the list of pivots on this layer
            std::pair<float, unsigned int> topOfQueue;
            while (copy_activePivots[layerIndex].size() > 0) {
                topOfQueue = copy_activePivots[layerIndex].top();
                copy_activePivots[layerIndex].pop();

                // get pivot information
                unsigned int const pivotIndex = topOfQueue.second;
                float const distance_Qp = topOfQueue.first;
                float const radius = pivotLayers[layerIndex].get_maxChildDistance(pivotIndex);

                // check if domain may contain a closer point
                if (distance_Qp <= dmin + radius) {
                    //
                    // iterate through the pivot domain
                    std::vector<unsigned int> const &pivotDomain =
                        pivotLayers[layerIndex].get_pivotChildren(pivotIndex);
                    for (int it1 = 0; it1 < pivotDomain.size(); it1++) {
                        unsigned int const childIndex = pivotDomain[it1];
                        float const distance_Qc = getDistance(queryIndex, childIndex, sparseMatrix);

                        // check if this child is the new closest point
                        if (distance_Qc < dmin) {
                            if (std::find(neighbors.begin(), neighbors.end(), childIndex) == neighbors.end()) {
                                dmin = distance_Qc;
                                index1 = childIndex;
                            }
                        }

                        // check if domain of child may contain a closer point
                        if (layerIndex < numLayers - 2) {
                            float const childRadius = pivotLayers[layerIndex + 1].get_maxChildDistance(childIndex);
                            if (distance_Qc <= dmin + childRadius) {
                                // add the pivot to the list
                                copy_activePivots[layerIndex + 1].push(std::make_pair(distance_Qc, childIndex));
                            }
                        }
                    }
                }
            }
        }

        //> Check the Domains of Intermediate Pivots Top-Down
        //  - improvement to make: if we have to examine the entire domain, then remove parent pivot from list
        //  - improvement to make: we can check if the entire domain is safe/invalid, then move to active list
        std::vector<PriorityQueue> copy_intermediatePivots = intermediatePivots;
        for (int layerIndex = 0; layerIndex < numLayers - 1; layerIndex++) {
            // printf("- IPD: Layer %d...\n",layerIndex);
            //
            // iterate through the list of pivots on this layer
            std::pair<float, unsigned int> topOfQueue;
            while (copy_intermediatePivots[layerIndex].size() > 0) {
                topOfQueue = copy_intermediatePivots[layerIndex].top();
                copy_intermediatePivots[layerIndex].pop();

                // get pivot information
                unsigned int const pivotIndex = topOfQueue.second;
                float const distance_Qp = topOfQueue.first;
                float const radius = pivotLayers[layerIndex].get_maxChildDistance(pivotIndex);

                // check if domain may contain a closer point
                if (distance_Qp <= dmin + radius) {
                    //
                    // iterate through the pivot domain
                    std::vector<unsigned int> const &pivotDomain =
                        pivotLayers[layerIndex].get_pivotChildren(pivotIndex);
                    for (int it1 = 0; it1 < pivotDomain.size(); it1++) {
                        unsigned int const childIndex = pivotDomain[it1];
                        float const distance_Qc = getDistance(queryIndex, childIndex, sparseMatrix);

                        // check if this child is the new closest point
                        if (distance_Qc < dmin) {
                            if (std::find(neighbors.begin(), neighbors.end(), childIndex) == neighbors.end()) {

                                // pivots are intermediate: must check that they are not invalidated
                                bool flag_validated = validatePoint(childIndex, queryIndex, neighbors, sparseMatrix);
                                if (flag_validated) {
                                    dmin = distance_Qc;
                                    index1 = childIndex;
                                }
                            }
                        }

                        // add pivot to examine if close enough
                        if (layerIndex < numLayers - 2) {
                            float const childRadius = pivotLayers[layerIndex + 1].get_maxChildDistance(childIndex);
                            if (distance_Qc <= dmin + childRadius) {
                                copy_intermediatePivots[layerIndex + 1].push(std::make_pair(distance_Qc, childIndex));
                            }
                        }
                    }
                }
            }
        }
        // printf("- Done searching...\n");


        if (dmin > 10000) break;                // large, magic number signaling no new neighbor possible
        unsigned int const x1 = index1;
        float const distance_Q1 = getDistance(queryIndex, x1, sparseMatrix);
        neighbors.push_back(x1);
        sparseMatrix._addNewReference(x1);  // storing distances from x1 to all other points



        /**
         * =================================================================
         *   STEP 2: Validation of Active Points
         *     - remove points that are invalidated by the new HSP neighbor
         *
         * =================================================================
         */
        
        // printf("-  Begin Validation...\n");
        
        //> Iterate through the intermediate list in each layer, top-down. 
        for (int layerIndex = 0; layerIndex < numLayers; layerIndex++) {
            // printf("- IPV: Layer %d...\n",layerIndex);

            // start by iterating through the intermediate pivots
            PriorityQueue list_copy = intermediatePivots[layerIndex];
            intermediatePivots[layerIndex] = PriorityQueue();
            std::pair<float, unsigned int> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                list_copy.pop();

                // get pivot information
                unsigned int const pivotIndex = topOfQueue.second;
                float const distance_Q2 = topOfQueue.first;
                float radius = 0;
                if (layerIndex < numLayers - 1) {
                    radius = pivotLayers[layerIndex].get_maxChildDistance(pivotIndex);
                }

                //> Pivot sufficiently close to Q. Add domain to lower layer. 
                if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {

                    //  - add the entire domain to the intermeidate list of the lower layer
                    if (layerIndex < numLayers - 1) {
                        std::vector<unsigned int> const &pivotDomain =
                            pivotLayers[layerIndex].get_pivotChildren(pivotIndex);
                        for (int it1 = 0; it1 < pivotDomain.size(); it1++) {
                            unsigned int const childIndex = pivotDomain[it1];
                            float const distance_Qc = getDistance(queryIndex, childIndex, sparseMatrix);
                            intermediatePivots[layerIndex+1].push(std::make_pair(distance_Qc, childIndex));
                        }
                    } 

                    //  - add the pivot back to the list
                    else {
                        intermediatePivots[layerIndex].push(std::make_pair(distance_Q2, pivotIndex));
                    }
                }

                //> Check pivot for invalidation
                else {
                    float const distance_12 = getDistance(x1, pivotIndex, sparseMatrix);

                    //  - proposition: entire domain invalidated
                    if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                        continue;
                    }
                    //  - add the pivot back to the list
                    else {
                        intermediatePivots[layerIndex].push(std::make_pair(distance_Q2, pivotIndex));
                    }
                }
            }
        }

         //> Iterate through the intermediate list in each layer, top-down. 
        for (int layerIndex = 0; layerIndex < numLayers; layerIndex++) {
            // printf("- APV: Layer %d...\n",layerIndex);

            // start by iterating through the intermediate pivots
            PriorityQueue list_copy = activePivots[layerIndex];
            activePivots[layerIndex] = PriorityQueue();
            std::pair<float, unsigned int> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                list_copy.pop();

                // get pivot information
                unsigned int const pivotIndex = topOfQueue.second;
                float const distance_Q2 = topOfQueue.first;
                float radius = 0;
                if (layerIndex < numLayers - 1) {
                    radius = pivotLayers[layerIndex].get_maxChildDistance(pivotIndex);
                }

                //> Pivot sufficiently close to Q. Add domain to lower layer list. 
                if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {

                    //  - add the entire domain to the intermeidate list of the lower layer
                    if (layerIndex < numLayers - 1) {
                        std::vector<unsigned int> const &pivotDomain =
                            pivotLayers[layerIndex].get_pivotChildren(pivotIndex);
                        for (int it1 = 0; it1 < pivotDomain.size(); it1++) {
                            unsigned int const childIndex = pivotDomain[it1];
                            float const distance_Qc = getDistance(queryIndex, childIndex, sparseMatrix);
                            intermediatePivots[layerIndex+1].push(std::make_pair(distance_Qc, childIndex));
                        }
                    } 

                    //  - add the pivot back to the list
                    else {
                        activePivots[layerIndex].push(std::make_pair(distance_Q2, pivotIndex));
                    }
                }

                //> Check pivot for invalidation
                else {
                    float const distance_12 = getDistance(x1, pivotIndex, sparseMatrix);

                    //  - proposition: entire domain invalidated
                    if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                        continue;
                    }
                    //  - proposition: Entire Domain Safe
                    else if (distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                        activePivots[layerIndex].push(std::make_pair(distance_Q2, pivotIndex));
                    }
                    //  - add the pivot back to the list
                    else {
                        intermediatePivots[layerIndex].push(std::make_pair(distance_Q2, pivotIndex));
                    }
                }
            }
        }

        //  - count the total number of points in lists
        unsigned int remainingPoints = 0;
        for (int layerIndex = 0; layerIndex < numLayers; layerIndex++) {
            remainingPoints += (activePivots[layerIndex].size() + intermediatePivots[layerIndex].size());
        }
    }

    return;
}