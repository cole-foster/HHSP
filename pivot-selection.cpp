#include "pivot-selection.hpp"

/**
 * @brief Index Construction: 
 * 
 * @param radius 
 * @param sparseMatrix 
 * @param pivotLayer 
 */
void PivotSelection::selectPivots(float radius, SparseMatrix &sparseMatrix, PivotLayer &pivotLayer) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    pivotLayer = PivotLayer(0, 2, radius);

    // incrementally add each point to the hierarchy
    for (unsigned int queryIndex = 0; queryIndex < datasetSize; queryIndex++) {
        float closestParentDistance = HUGE_VAL;
        unsigned int closestParent;

        // compute distance to all existing parents
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = pivotLayer.pivotIndices->begin(); it1 != pivotLayer.pivotIndices->end(); it1++) {
            unsigned int const pivotIndex = (*it1);
            float const distance = sparseMatrix._computeDistance(queryIndex,pivotIndex);
            if (distance <= radius) {
                if (distance < closestParentDistance) {
                    closestParentDistance = distance;
                    closestParent = pivotIndex;
                }
            }
        }

        // add query to domain of the closest parent 
        if (closestParentDistance <= 10000) { // some magically large number
            pivotLayer.addPivotChild(closestParent,queryIndex);
        } 

        // if no parent, query becomes a pivot
        else {
            pivotLayer.addPivot(queryIndex);
            pivotLayer.addPivotChild(queryIndex,queryIndex);
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
bool PivotSelection::validatePivotSelection(PivotLayer &pivotLayer, SparseMatrix &sparseMatrix) {
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
        tsl::sparse_set<unsigned int> const& pivotDomain = pivotLayer.get_pivotChildren(pivotIndex);

        tsl::sparse_set<unsigned int>::const_iterator it2;
        for (it2 = pivotDomain.begin(); it2 != pivotDomain.end(); it2++) {
            unsigned int const childIndex = (*it2);

            // check within radius
            float const distance = sparseMatrix._computeDistance(pivotIndex,childIndex);
            if (distance > radius) {
                printf("Point not within radius of pivot!\n");
                printf("p=%u,c=%u,d(p,c)=%.4f\n",pivotIndex,childIndex,distance);
                return false;
            }

            // remove child from list, ensure not seen before
            std::vector<unsigned int>::iterator it3 = std::find(pointList.begin(),pointList.end(),childIndex);
            if (it3 == pointList.end()) {
                printf("Child has already been represented! Coverage not minimal!\n");
                printf("c=%u\n",childIndex);
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

void PivotSelection::printSet(tsl::sparse_set<unsigned int> const &set) {
    printf("{");
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}
