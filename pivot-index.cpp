// Cole Foster
// 2023-03-20
#include "pivot-index.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

/**
 * @brief Index Construction:
 *
 * @param radius
 * @param sparseMatrix
 * @param pivotLayer
 */
void PivotIndex::Greedy2L(float radius, SparseMatrix& sparseMatrix, std::vector<Pivot>& pivot_list) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    pivot_list.clear();

    // incrementally add each point to the hierarchy
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        float closestParentDistance = HUGE_VAL;
        Pivot* closestParent;

        // compute distance to all existing parents
        std::vector<Pivot>::iterator it1;
        for (it1 = pivot_list.begin(); it1 != pivot_list.end(); it1++) {
            Pivot& pivot = (*it1);
            float const distance = sparseMatrix._computeDistance(queryIndex, pivot._index);
            if (distance <= radius) {
                if (distance < closestParentDistance) {
                    closestParentDistance = distance;
                    closestParent = &pivot;
                }
            }
        }

        // add query to domain of the closest parent
        if (closestParentDistance <= 10000) {                            // some magically large number
            closestParent->addChild(queryIndex, closestParentDistance);  // added as a point
        }

        // if no parent, query becomes a pivot
        else {
            Pivot queryPivot = Pivot(queryIndex, radius);
            queryPivot.addChild(queryIndex, 0.0f);  // added as a point
            pivot_list.push_back(queryPivot);
        }
    }

    return;
}

/**
 * @brief Construct the Cover Tree incrementally.
 *
 * Concepts:
 * - A query Q has a parent p in layer ell if d(Q,p) < r^{ell} - r^{ell+1}
 * - A pivot p in layer ell may recursively have children that are parents of Q in lower layers if
 *   d(Q,p) < r^{ell}
 * - We need to find the lowest layer without any parents of Q (otherwise this invalidates the separation invariant)
 */
void PivotIndex::Greedy_MultiLayer(std::vector<float> radiusVector, SparseMatrix& sparseMatrix,
                                   std::vector<Pivot>& pivot_list, std::vector<unsigned int>& pivotsPerLayer) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    radiusVector.push_back(0);
    int numberOfLayers = radiusVector.size();
    pivotsPerLayer.resize(numberOfLayers-1);
    pivot_list.clear();

    // incremental construction in a top-down fashion
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        int orphanLayer = 0;  // 1 below lowest layer with a parent
        Pivot* lowestParent; // lowest parent, orphanLayer-1

        // top-down, find lowest layer with a parent
        std::vector<Pivot*> grandParents{};
        for (int layerIndex = 0; layerIndex < numberOfLayers - 1; layerIndex++) {

            //> Collect the Spotlight: nodes in layer local to Q
            std::vector<Pivot*> spotlight{}; 
            if (layerIndex == 0) { // collect pointers to all top layer pivots
                for (int it1 = 0; it1 < (int) pivot_list.size(); it1++) {
                    spotlight.push_back(&pivot_list[it1]);
                }
            } else { // collect pointers to children of grandparents
                for (int it1 = 0; it1 < (int) grandParents.size(); it1++) {
                    Pivot* grandPivot = grandParents[it1];
                    for (int it2 = 0; it2 < (int) grandPivot->_pivotDomain.size(); it2++) {
                        spotlight.push_back(&grandPivot->_pivotDomain[it2]);
                    }
                }
            }

            //> Find the closest parent of Q from Spotlight
            float closestParentDistance = HUGE_VAL;
            grandParents.clear();
            for (int it2 = 0; it2 < (int) spotlight.size(); it2++) {
                Pivot* candidateParentPivot = spotlight[it2];
                float const distance = sparseMatrix._computeDistance(queryIndex, candidateParentPivot->_index);

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

        //> Add query to the tree
        Pivot queryPivot = Pivot(queryIndex, radiusVector[orphanLayer], orphanLayer);
        Pivot* queryPivotPtr;
        if (orphanLayer == 0) {
            pivot_list.push_back(queryPivot);
            queryPivotPtr = &pivot_list.back();
        } else {
            lowestParent->addChild(queryPivot, sparseMatrix._computeDistance(queryIndex,lowestParent->_index));
            queryPivotPtr = &lowestParent->_pivotDomain.back();
        }
        
        //> Add query as a child of itself in all lower layers
        for (int layerIndex = orphanLayer; layerIndex < numberOfLayers-1; layerIndex++) {
            queryPivot = Pivot(queryIndex, radiusVector[layerIndex+1], layerIndex+1);
            queryPivotPtr->addChild(queryPivot,0.0f);
            queryPivotPtr = &queryPivotPtr->_pivotDomain.back();
            pivotsPerLayer[layerIndex]++;
        }
    }

    return;
}






bool recursivePivotDomainValidation(Pivot const& pivot, int const numberOfLayers, SparseMatrix& sparseMatrix, std::vector<unsigned int>& pointList) {

    // iterate through pivot domain, check if valid
    float const radius = pivot._radius;
    std::vector<Pivot> const& pivotDomain = pivot._pivotDomain;
    std::vector<Pivot>::const_iterator it2;
    for (it2 = pivotDomain.begin(); it2 != pivotDomain.end(); it2++) {
        Pivot const& child = (*it2);

        // check within radius
        float const distance = sparseMatrix._computeDistance(pivot._index, child._index);
        if (distance > radius) {
            printf("Point not within radius of pivot!\n");
            printf("p=%u,c=%u,d(p,c)=%.4f\n", pivot._index, child._index, distance);
            return false;
        }

        // if bottom layer, remove from list
        if (child._level == numberOfLayers-1) {
            // remove child from list, ensure not seen before
            std::vector<unsigned int>::iterator it3 = std::find(pointList.begin(), pointList.end(), child._index);
            if (it3 == pointList.end()) {
                printf("Child has already been represented! Coverage not minimal!\n");
                printf("c=%u\n", child._index);
                return false;
            } else {
                pointList.erase(it3);
            }
        } 

        // otherwise, check their domains
        else {
            if (!recursivePivotDomainValidation(child, numberOfLayers, sparseMatrix, pointList)) {
                return false;
            }
        }
    }

    return true;
}




/**
 * @brief ensuring minimal coverage of the set
 *
 * @param sparseMatrix
 * @param pivotLayerVector
 */
bool PivotIndex::validatePivotSelection(std::vector<Pivot>& pivotsList, int const numberOfLayers, SparseMatrix& sparseMatrix) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;

    std::vector<unsigned int> pointList{};
    pointList.resize(datasetSize);
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        pointList[queryIndex] = queryIndex;
    }

    // now, iterate through the pivot domains and remove points one-by-one
    std::vector<Pivot>::const_iterator it1;
    for (it1 = pivotsList.begin(); it1 != pivotsList.end(); it1++) {
        Pivot const& pivot = (*it1);
        if (!recursivePivotDomainValidation(pivot, numberOfLayers, sparseMatrix, pointList)) {
            return false;
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


void PivotIndex::printSet(std::vector<unsigned int> const &set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}