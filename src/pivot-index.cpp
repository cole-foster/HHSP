#include "pivot-index.hpp"

/**
 * @brief Greedy Construction of the Cover Tree with Finite Layers (Defined by Radii). Greedy as in built incrementally
 * and does not perform nearest-parent assignment.
 *
 * @param radiusVector
 * @param sparseMatrix
 * @param coverTree
 */
void PivotIndex::CoverTree_Greedy(std::vector<float> radiusVector, SparseMatrix& sparseMatrix,
                                  std::vector<PivotLayer>& coverTree) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    radiusVector.push_back(0);
    int const numLayers = radiusVector.size();

    //> Change to Effective Radius Vector -> Makes them Covering Radii for entire pivot domains
    std::vector<float> temp_vec = radiusVector;
    for (int i = 0; i < (int)radiusVector.size(); i++) {
        for (int j = i + 1; j < (int)temp_vec.size(); j++) {
            radiusVector[i] += (float)temp_vec[j];
        }
    }

    //> Initialize Each Layer of the Cover Tree
    coverTree.clear();
    coverTree.resize(numLayers - 1);
    for (int layerIndex = 0; layerIndex < numLayers - 1; layerIndex++) {
        coverTree[layerIndex] = PivotLayer(layerIndex, numLayers, radiusVector[layerIndex]);
    }

    //> Incremental construction in a top-down fashion
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        // printf("Q: %u...\n",queryIndex);
        int orphanLayer = 0;        // 1 Below Lowest Layer With a Parent
        unsigned int lowestParent;  // Closest Parent in the Layer Above Orphanage

        //> Top-Down approach, find the lowest layer with a parent
        std::vector<unsigned int> grandParents{};
        for (int layerIndex = 0; layerIndex <= numLayers - 1; layerIndex++) {
            //> Collect the Spotlight |--> Pivots in layer local to Q (may be a parent)
            std::vector<unsigned int> spotlight{};
            if (layerIndex == 0) {
                // Collect All Top Layer Pivots
                tsl::sparse_set<unsigned int> const& topLayerPivots = *(coverTree[0].get_pivotIndices_ptr());
                spotlight.insert(spotlight.end(), topLayerPivots.begin(), topLayerPivots.end());

            } else {
                // Collect the Children of "Grandparents", may be parents of Q
                for (int it1 = 0; it1 < (int)grandParents.size(); it1++) {
                    unsigned int const grandPivot = grandParents[it1];
                    std::vector<unsigned int> const& pivotDomain =
                        coverTree[layerIndex - 1].get_pivotChildren(grandPivot);
                    spotlight.insert(spotlight.end(), pivotDomain.begin(), pivotDomain.end());
                }
            }

            //> Find the closest parent of Q from the Spotlight
            float closestParentDistance = HUGE_VAL;
            grandParents.clear();  // grandparents: any pivot within the covering radius r of Q
            for (int it2 = 0; it2 < (int)spotlight.size(); it2++) {
                unsigned int const candidateParentPivot = spotlight[it2];
                float const distance = sparseMatrix._computeDistance(queryIndex, candidateParentPivot);

                // Test if Q is within the Covering Radius
                if (distance <= radiusVector[layerIndex]) {
                    grandParents.push_back(candidateParentPivot);  // could have children as lower parent

                    // Test if Q is an actual parent of Q
                    if (distance <= radiusVector[layerIndex] - radiusVector[layerIndex + 1]) {
                        // is this new lowest parent layer?
                        if (orphanLayer < layerIndex + 1) {
                            orphanLayer = layerIndex + 1;  // no parent yet in layer below
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
            if (grandParents.empty()) break;
        }

        //> Nesting: Add Query To Orphan Layer and All Layers Below
        for (int layerIndex = orphanLayer; layerIndex < numLayers; layerIndex++) {
            //
            //> Becomes a Pivot on Every Lower Layer Except Bottom: Not Explicitly Represented
            if (layerIndex < numLayers - 1) {
                coverTree[layerIndex].addPivot(queryIndex);
            }

            //> Becomes a Child Of Itself
            if (layerIndex > 0) {
                coverTree[layerIndex - 1].addChild(queryIndex, queryIndex);
            }

            //> Becomes a Parent Of Itself
            if (layerIndex < numLayers - 1 && layerIndex > orphanLayer) {
                coverTree[layerIndex].addParent(queryIndex, queryIndex);
            }
        }

        // Add The Parent of the Query, Recursively Update Parents MaxChildDistance
        if (orphanLayer > 0) {
            if (orphanLayer < numLayers - 1) {
                coverTree[orphanLayer].addParent(queryIndex, lowestParent);
            }
            coverTree[orphanLayer - 1].addChild(lowestParent, queryIndex);
            coverTree[orphanLayer - 1].updateMaxChildDistance(lowestParent,
                                                              sparseMatrix._computeDistance(lowestParent, queryIndex));

            // Update Max Child Distance for All Parents of Parents
            unsigned int parent = lowestParent;
            for (int layerIndex = orphanLayer - 1; layerIndex > 0; layerIndex--) {
                unsigned int grandParent = coverTree[layerIndex].get_parentIndex(parent);
                coverTree[layerIndex - 1].updateMaxChildDistance(
                    grandParent, sparseMatrix._computeDistance(grandParent, queryIndex));
                parent = grandParent;
            }
        }
    }

    return;
}

// used for cover tree validation, below
bool recursiveDepthFirstCheck(std::vector<unsigned int> parents, int const layerIndex,
                              std::vector<PivotLayer>& coverTree, SparseMatrix& sparseMatrix,
                              std::vector<unsigned int>& pointList);

/**
 * @brief Validate the Coverage and Nesting Properties of the Cover Tree, but not Separation.
 *
 * @param coverTree
 * @param sparseMatrix
 * @return true
 * @return false
 */
bool PivotIndex::validateCoverTree(std::vector<PivotLayer>& coverTree, SparseMatrix& sparseMatrix) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    int const numLayers = coverTree.size() + 1;

    // print stats:
    printf("LayerID, Radius, NumPivots\n");
    for (int layerIndex = 0; layerIndex < numLayers - 1; layerIndex++) {
        printf("%d, %.4f, %u\n", layerIndex, coverTree[layerIndex].radius, coverTree[layerIndex].pivotIndices->size());
    }
    printf("%d, %.4f, %u\n", numLayers - 1, 0.0000, datasetSize);

    // each point should be represented
    std::vector<unsigned int> pointList{};
    pointList.resize(datasetSize);
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        pointList[queryIndex] = queryIndex;
    }

    // Top-Layer Iteration
    tsl::sparse_set<unsigned int> const& pivotIndices = *(coverTree[0].get_pivotIndices_ptr());
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivotIndex = (*it1);
        std::vector<unsigned int> parentsList{};
        parentsList.resize(numLayers);
        parentsList[0] = pivotIndex;
        bool flag_success = recursiveDepthFirstCheck(parentsList, 0, coverTree, sparseMatrix, pointList);
        if (!flag_success) return false;
    }

    // should be no points left
    if (pointList.size() > 0) {
        printf("There are points not covered by the index! Coverage violated!\n");
        printf("  -> %u: ", pointList.size());
        PivotIndex::printSet(pointList);
        return false;
    }

    return true;
}

/**
 * @brief Recursive Function for Performing Depth-First Cover Tree Validation
 *
 * @param parents
 * @param layerIndex
 * @param coverTree
 * @param sparseMatrix
 * @param pointList
 * @return true
 * @return false
 */
bool recursiveDepthFirstCheck(std::vector<unsigned int> parents, int const layerIndex,
                              std::vector<PivotLayer>& coverTree, SparseMatrix& sparseMatrix,
                              std::vector<unsigned int>& pointList) {
    int numLayers = coverTree.size() + 1;

    // Iterate through the Domain of Lowest Parent
    std::vector<unsigned int> const& pivotDomain = coverTree[layerIndex].get_pivotChildren(parents[layerIndex]);
    for (int it1 = 0; it1 < (int)pivotDomain.size(); it1++) {
        unsigned int const childIndex = pivotDomain[it1];
        float childRadius = 0;       // radius of the child pivot
        float childMaxDistance = 0;  // max distance to child of child
        if (layerIndex + 1 < numLayers - 1) {
            childRadius = coverTree[layerIndex + 1].radius;
            // childMaxDistance = coverTree[layerIndex + 1].get_maxChildDistance(childIndex);
        }

        // check childIndex within radius/maxDistance of every parent index
        for (int coarserLayerIndex = 0; coarserLayerIndex <= layerIndex; coarserLayerIndex++) {
            unsigned int const coarserParent = parents[coarserLayerIndex];
            float const pivotRadius = coverTree[coarserLayerIndex].radius;
            float const pivotMaxDistance = coverTree[coarserLayerIndex].get_maxChildDistance(coarserParent);

            float const distance = sparseMatrix._computeDistance(coarserParent, childIndex);
            if (distance > pivotRadius - childRadius) {
                printf("ChildIndex %u on L-%d not within radius of ParentIndex %u on L-%d!\n", childIndex,
                       layerIndex + 1, coarserParent, coarserLayerIndex);
                return false;
            }

            if (distance > pivotMaxDistance) {
                printf(
                    "ChildIndex %u on L-%d -max child distance- not within max-child-distance of ParentIndex %u on "
                    "L-%d!\n",
                    childIndex, layerIndex + 1, coarserParent, coarserLayerIndex);
                return false;
            }
        }

        if (layerIndex + 1 < numLayers - 1) {
            parents[layerIndex + 1] = childIndex;
            recursiveDepthFirstCheck(parents, layerIndex + 1, coverTree, sparseMatrix, pointList);
        } else {
            // if bottom layer, remove points from list
            std::vector<unsigned int>::iterator it3a = std::find(pointList.begin(), pointList.end(), childIndex);
            if (it3a == pointList.end()) {
                printf("Child has already been represented! Coverage not minimal!\n");
                printf("c=%u\n", childIndex);
                return false;
            } else {
                pointList.erase(it3a);
            }
        }
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
