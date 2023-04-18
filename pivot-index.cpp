#include "pivot-index.hpp"


void PivotIndex::Greedy_3L(std::vector<float> radiusVector, SparseMatrix& sparseMatrix, std::vector<PivotLayer>& pivotLayers) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    if (radiusVector.size() != 2) {
        printf("Incorrect Radius Vector for 3L\n");
        return;
    }
    radiusVector.push_back(0);

    //> Initialize Pivot Hierarchy
    pivotLayers.clear();
    pivotLayers.resize(2);
    pivotLayers[0] = PivotLayer(0,3,radiusVector[0]);
    pivotLayers[1] = PivotLayer(1,3,radiusVector[1]);

    //> Incremental construction in a top-down fashion
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        int orphanLayer = 0;  // 1 Below Layer without Parent
        unsigned int lowestParent;  // Closest Parent in layer above orphanage

        // top-down, find lowest layer with a parent
        std::vector<unsigned int> grandParents{};
        for (int layerIndex = 0; layerIndex <= 1; layerIndex++) {

            //> Collect the Spotlight |--> Pivots in layer local to Q
            std::vector<unsigned int> spotlight{}; 
            if (layerIndex == 0) { // collect all top layer pivots
                tsl::sparse_set<unsigned int> const& topLayerPivots = *(pivotLayers[0].get_pivotIndices_ptr());
                spotlight.insert(spotlight.end(), topLayerPivots.begin(), topLayerPivots.end());

            } else { // collect children of grandparents
                for (int it1 = 0; it1 < (int) grandParents.size(); it1++) {
                    unsigned int const grandPivot = grandParents[it1];
                    std::vector<unsigned int> const& pivotDomain = pivotLayers[layerIndex-1].get_pivotChildren(grandPivot);
                    spotlight.insert(spotlight.end(), pivotDomain.begin(), pivotDomain.end());
                }
            }

            //> Find the closest parent of Q from Spotlight
            float closestParentDistance = HUGE_VAL;
            grandParents.clear();
            for (int it2 = 0; it2 < (int) spotlight.size(); it2++) {
                unsigned int const candidateParentPivot = spotlight[it2];
                float const distance = sparseMatrix._computeDistance(queryIndex, candidateParentPivot);

                // test if pivot may cover a future parent
                if (distance <= radiusVector[layerIndex]) {
                    grandParents.push_back(candidateParentPivot);  // could have children as lower parent

                    // test if pivot is indeed a parent
                    if (distance <= radiusVector[layerIndex] - radiusVector[layerIndex + 1]) {

                        // is this new lowest parent layer?
                        if (orphanLayer < layerIndex + 1) {
                            orphanLayer = layerIndex + 1;
                        }

                        // is this the closest parent?
                        if (distance < closestParentDistance) {
                            closestParentDistance = distance;
                            lowestParent = candidateParentPivot;
                        }
                    }
                }
            }

            //> Break if no potential future parents
            if (grandParents.empty()) {
                break;
            }
        }
        
        // If query a top-layer pivot
        if (orphanLayer == 0) {

            // add query to top layer
            pivotLayers[0].addPivot(queryIndex);
            
            // add query to second layer
            pivotLayers[1].addPivot(queryIndex);
            pivotLayers[1].addParent(queryIndex, queryIndex);
            pivotLayers[0].addChild(queryIndex, queryIndex);
            pivotLayers[0].updateMaxChildDistance(queryIndex, 0.0f);

            // add query to bottom layer
            pivotLayers[1].addChild(queryIndex, queryIndex);
            pivotLayers[1].updateMaxChildDistance(queryIndex, 0.0f);
            pivotLayers[0].updateMaxChildDistance(queryIndex, 0.0f);
        }

        // If query a second layer pivot
        else if (orphanLayer == 1) {

            // add query to second layer
            pivotLayers[1].addPivot(queryIndex);
            pivotLayers[1].addParent(queryIndex, lowestParent);
            pivotLayers[0].addChild(lowestParent, queryIndex);
            pivotLayers[0].updateMaxChildDistance(lowestParent, sparseMatrix._computeDistance(lowestParent, queryIndex));

            // add query to bottom layer
            pivotLayers[1].addChild(queryIndex, queryIndex);
            pivotLayers[1].updateMaxChildDistance(queryIndex, 0.0f);
            pivotLayers[0].updateMaxChildDistance(lowestParent, 0.0f);
        }


        // If query a bottom layer pivot only
        else {

            // add query to bottom layer
            pivotLayers[1].addChild(lowestParent, queryIndex);
            pivotLayers[1].updateMaxChildDistance(lowestParent, sparseMatrix._computeDistance(lowestParent, queryIndex));
            unsigned int const grandParent = pivotLayers[1].get_parentIndex(lowestParent);
            pivotLayers[0].updateMaxChildDistance(grandParent, sparseMatrix._computeDistance(grandParent, queryIndex));
        }

    }

    return;
}


/**
 * @brief Construct a 2-Layer Cover Tree
 *
 * @param radius
 * @param sparseMatrix
 * @param pivotLayer
 */
void PivotIndex::Greedy_2L(float radius, SparseMatrix& sparseMatrix, std::vector<PivotLayer>& pivotLayers) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    pivotLayers.clear();
    pivotLayers.push_back(PivotLayer(0, 2, radius));

    // incrementally add each point to the hierarchy
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        float closestParentDistance = HUGE_VAL;
        unsigned int closestParent;

        // compute distance to all existing parents
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = pivotLayers[0].pivotIndices->begin(); it1 != pivotLayers[0].pivotIndices->end(); it1++) {
            unsigned int const pivotIndex = (*it1);
            float const distance = sparseMatrix._computeDistance(queryIndex, pivotIndex);
            if (distance <= radius) {
                if (distance < closestParentDistance) {
                    closestParentDistance = distance;
                    closestParent = pivotIndex;
                }
            }
        }

        // add query to domain of the closest parent
        if (closestParentDistance <= 10000) {  // some magically large number
            pivotLayers[0].addChild(closestParent, queryIndex);
            pivotLayers[0].updateMaxChildDistance(closestParent, closestParentDistance);
        }

        // if no parent, query becomes a pivot
        else {
            pivotLayers[0].addPivot(queryIndex);
            pivotLayers[0].addChild(queryIndex, queryIndex);
        }
    }

    return;
}





/**
 * @brief ensuring minimal coverage of the set
 *
 * @param sparseMatrix
 * @param pivotLayerVector
 */
bool PivotIndex::validatePivotSelection(PivotLayer& pivotLayer, SparseMatrix& sparseMatrix) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    float const radius = pivotLayer.radius;
    std::vector<unsigned int> pointList{};
    pointList.resize(datasetSize);
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        pointList[queryIndex] = queryIndex;
    }

    // now, iterate through the pivot domains and remove points one-by-one
    tsl::sparse_set<unsigned int> const& pivotIndices = *pivotLayer.get_pivotIndices_ptr();
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        std::vector<unsigned int> const& pivotDomain = pivotLayer.get_pivotChildren(pivotIndex);

        std::vector<unsigned int>::const_iterator it2;
        for (it2 = pivotDomain.begin(); it2 != pivotDomain.end(); it2++) {
            unsigned int const childIndex = (*it2);

            // check within radius
            float const distance = sparseMatrix._computeDistance(pivotIndex, childIndex);
            if (distance > radius) {
                printf("Point not within radius of pivot!\n");
                printf("p=%u,c=%u,d(p,c)=%.4f\n", pivotIndex, childIndex, distance);
                return false;
            }

            // remove child from list, ensure not seen before
            std::vector<unsigned int>::iterator it3 = std::find(pointList.begin(), pointList.end(), childIndex);
            if (it3 == pointList.end()) {
                printf("Child has already been represented! Coverage not minimal!\n");
                printf("c=%u\n", childIndex);
                return false;
            } else {
                pointList.erase(it3);
            }
        }
    }

    // should be no points left
    if (pointList.size() > 0) {
        printf("There are points not covered by the index! Coverage violated!\n");
        return false;
    }

    return true;
}

bool PivotIndex::validatePivotSelection_3L(std::vector<PivotLayer>& pivotLayers, SparseMatrix& sparseMatrix) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;

    // each point should be represented
    std::vector<unsigned int> pointList{};
    pointList.resize(datasetSize);
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        pointList[queryIndex] = queryIndex;
    }

    //> Depth-First Traversal, remove each point from list
    
    // Top-Layer Iteration
    tsl::sparse_set<unsigned int> const& pivotIndices = *(pivotLayers[0].get_pivotIndices_ptr());
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex1 = (*it1);
        float const radius1 = pivotLayers[0].radius;
        float const maxDist1 = pivotLayers[0].get_maxChildDistance(pivotIndex1);
        std::vector<unsigned int> const& pivotDomain1 = pivotLayers[0].get_pivotChildren(pivotIndex1);

        // Second-Layer Iteration
        std::vector<unsigned int>::const_iterator it2;
        for (it2 = pivotDomain1.begin(); it2 != pivotDomain1.end(); it2++) {
            unsigned int const pivotIndex2 = (*it2);
            std::vector<unsigned int> const& pivotDomain2 = pivotLayers[1].get_pivotChildren(pivotIndex2);
            float const radius2 = pivotLayers[1].radius;
            float const maxDist2 = pivotLayers[1].get_maxChildDistance(pivotIndex2);

            // check pivot2 within radius of pivot1
            float const distance12 = sparseMatrix._computeDistance(pivotIndex1,pivotIndex2);
            if (distance12 > radius1 || distance12 > maxDist1) {
                printf("L2 Pivot not within radius of L1 Pivot!\n");
                printf("p=%u,c=%u,d(p,c)=%.4f,r=%.4f,MCD=%.4f\n", pivotIndex1, pivotIndex2, distance12, radius1, maxDist1);
                return false;
            }

            // Second-Layer Iteration
            std::vector<unsigned int>::const_iterator it3;
            for (it3 = pivotDomain2.begin(); it3 != pivotDomain2.end(); it3++) {
                unsigned int const pivotIndex3 = (*it3);

                // check pivot3 within radius of pivot2
                float const distance23 = sparseMatrix._computeDistance(pivotIndex2,pivotIndex3);
                if (distance23 > radius2 || distance23 > maxDist2) {
                    printf("L3 Pivot not within radius of L2 Pivot!\n");
                    printf("p=%u,c=%u,d(p,c)=%.4f,r=%.4f,MCD=%.4f\n", pivotIndex2, pivotIndex3, distance23, radius2, maxDist2);
                    return false;
                }

                // check pivot3 within radius of pivot1
                float const distance13 = sparseMatrix._computeDistance(pivotIndex1,pivotIndex3);
                if (distance13 > radius1 || distance13 > maxDist1) {
                    printf("L3 Pivot not within radius of L1 Pivot!\n");
                    printf("p=%u,c=%u,d(p,c)=%.4f,r=%.4f,MCD=%.4f\n", pivotIndex1, pivotIndex3, distance13, radius1, maxDist1);
                    return false;
                }

                // now remove pivot3 from representation list
                std::vector<unsigned int>::iterator it3a = std::find(pointList.begin(), pointList.end(), pivotIndex3);
                if (it3a == pointList.end()) {
                    printf("Child has already been represented! Coverage not minimal!\n");
                    printf("c=%u\n", pivotIndex3);
                    return false;
                } else {
                    pointList.erase(it3a);
                }   
            }
        }
    }

    // should be no points left
    if (pointList.size() > 0) {
        printf("There are points not covered by the index! Coverage violated!\n");
        printf("  -> %u: ",pointList.size()); PivotIndex::printSet(pointList);
        return false;
    }

    return true;
}

void PivotIndex::printSet(std::vector<unsigned int> const& set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}
