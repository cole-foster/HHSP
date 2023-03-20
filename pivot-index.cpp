// Cole Foster
// 2023-03-20
#include "pivot-index.hpp"
#include <cmath>
#include <cstdio>
#include <algorithm>

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
 * @brief ensuring minimal coverage of the set
 *
 * @param sparseMatrix
 * @param pivotLayerVector
 */
bool PivotIndex::validatePivotSelection(std::vector<Pivot>& pivotsList, SparseMatrix& sparseMatrix) {
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
    }

    // should be no points left
    if (pointList.size() > 0) {
        printf("There are points not covered by the index! Coverage violated!\n");
        return false;
    }

    return true;
}

