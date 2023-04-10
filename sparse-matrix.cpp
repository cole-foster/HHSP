/*
Copyright 2022, Brown University, Providence, RI.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose other than its incorporation into a
commercial product or service is hereby granted without fee, provided
that the above copyright notice appear in all copies and that both
that copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
// Cole Foster
// 01-24-2022
#include "sparse-matrix.hpp"
#include <cmath>

SparseMatrix::SparseMatrix(float* const& dataPointer, unsigned int const datasetSize, unsigned int const dimension)
    : _dataPointer(dataPointer), _datasetSize(datasetSize), _dimension(dimension) {
}

/**
 * @brief compute distance between two indices
 * 
 * @param index1 
 * @param index2 
 * @return float const 
 */
float const SparseMatrix::_computeDistance(unsigned int const index1, unsigned int const index2) {
    _distanceComputationCount++;

    float distance = 0;
    for (unsigned int d = 0; d < _dimension; d++) {
        float difference = (_dataPointer[index1 * _dimension + d]) - (_dataPointer[index2 * _dimension + d]);
        distance += difference * difference;
    }

    return std::sqrt(distance);
}

/**
 * @brief Special Storage of Distance Computations
 * 
 * @param index 
 */
void SparseMatrix::_addNewReference(unsigned int index) {
    std::vector<float> new_vec(_datasetSize,-1);
    _distanceMatrix.insert({index,new_vec});
}

/**
 * @brief compute distance between two indices
 * 
 * @param index1 
 * @param index2 
 * @return float const 
 */
float const SparseMatrix::_getDistance(unsigned int const index1, unsigned int const index2) {
    
    // make sure this vector is in the map. if not, just compute and leave. 
    // COMMENT THIS OUT: INEFFICIENT, ONLY CALL WHEN IT IS IN THE MAP
    std::unordered_map<unsigned int, std::vector<float>>::const_iterator it1;
    it1 = _distanceMatrix.find(index1);
    if (it1 == _distanceMatrix.end()) {
        return _computeDistance(index1,index2);
        printf("Not in map loser\n");
    }  

    // grab distance from matrix
    float distance = _distanceMatrix[index1][index2];
    if (distance < 0) {
        // compute and store distance
        distance = _computeDistance(index1,index2);
        _distanceMatrix[index1][index2] = distance;
    }

    return distance;
}