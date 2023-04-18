#ifndef SparseMatrix_h
#define SparseMatrix_h
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

#include <cstddef>
#include <tsl/sparse_map.h>


/**
 * @brief sparse matrix to hold all pairwise distance computations computed. Stored as lower triangular matrix. Uses a vector of hashmaps for speed.
 *
 */
class SparseMatrix {
  public:
    SparseMatrix(){};
    SparseMatrix(float *const &dataPointer, unsigned int const datasetSize, unsigned int const dimension);
    ~SparseMatrix(){};
    float const _computeDistance(unsigned int const index1, unsigned int const index2);
    float const _getDistance(unsigned int const index1, unsigned int const index2);
    void _addNewReference(unsigned int index);
    inline void _clear() {_distanceMatrix.clear();}

    float *_dataPointer = NULL;
    unsigned int _datasetSize = 0;
    unsigned int _dimension = 0;
    unsigned long long int _distanceComputationCount = 0;

  private:
    // distance computation array
    tsl::sparse_map<unsigned int, std::vector<float>> _distanceMatrix{};

    
};

#endif  // SparseMatrix_h