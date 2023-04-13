#ifndef datasets_hpp
#define datasets_hpp

/*
Copyright 2019, Brown University, Providence, RI.
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
// November 29th, 2021

#include <string>

namespace Datasets {

// generate dataset of [N,D] floating point numbers uniformly distributed from [-1,1]
void uniformDataset(float*& dataPointer, unsigned int dimension, unsigned int N, float b, unsigned int seed);

void clusterDataset(float*& dataPointer, unsigned int dimension, unsigned int& datasetSize, std::string data_directory);

void datasetShuffle(float*& dataPointer, unsigned int dimension, unsigned int datasetSize);

void extractTestset(unsigned int dimension, float*& dataPointer, unsigned int& datasetSize, float*& testPointer, unsigned int testsetSize);

void Cities(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory);

void LA(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, 
           float*& testPointer, unsigned int& testSetSize, 
           std::string dataDirectory); 


};  // namespace Datasets



#endif  // datasets_hpp