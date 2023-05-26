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

// class to create random datasets
// uses random seed for reproductibility
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
#include <fstream>


namespace Datasets {
std::string data_directory = "/users/cfoste18/data/cfoste18/datasets/";

// generate dataset of [N,D] floating point numbers
// uniformly distributed from [-1,1]
void uniformDataset(float*& dataPointer, unsigned int dimension, unsigned int datasetSize, unsigned int seed) {
    srand(seed);  // for same initializations of random variables

    unsigned long long int numElements = (unsigned long long int)dimension * (unsigned long long int) datasetSize;
    dataPointer = new float[numElements];

    for (unsigned long long int i = 0; i < numElements; i++) {
        dataPointer[i] = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * 2 - 1;
    }
    
    return;
}

/**
 * @brief Get a clustered dataset
 *
 * @param datasetSize           desired datasetSize, adjusted value
 * @param dimension             dimension of vectors
 * @param numClusters           number of clusters in dataset
 * @param variance              variance of the gaussian perturbations
 * @return std::vector<std::vector<double>>
 */
void clusteredDataset(float*& dataPointer, unsigned int dimension, unsigned int &datasetSize, unsigned int seed, int const numClusters,
                      float const variance) {
    srand(seed);
    unsigned int clusterSize = (unsigned int) ceil((double)(datasetSize) / (double)numClusters);
    datasetSize = (clusterSize * (unsigned int) numClusters);

    // initialize the vectors with length this
    unsigned long long int numElements = (unsigned long long int) datasetSize * (unsigned long long int) dimension;
    dataPointer = new float[numElements];

    // shuffling the dataset for a random order
    std::mt19937 rng(11);
    std::vector<unsigned int> indices(datasetSize, 0);
    std::iota(indices.begin(), indices.end(), 0);       // 0 , 1, 2, ... datasetSize-1
    std::shuffle(indices.begin(), indices.end(), rng);  // shuffle indices by rng seed

    // prepping for the gassians
    std::default_random_engine generator;

    // Add each cluster
    unsigned long long int UL_DIMENSION = (unsigned long long int) dimension;
    unsigned long long int UL_CLUSTERSIZE = (unsigned long long int) clusterSize;
    std::vector<float> center_point; 
    center_point.resize(dimension);
    unsigned long long int c, d, i;
    for (c = 0; c < (unsigned long long int) numClusters; c++) {
        
        // get the random cluster center
        center_point.resize(dimension);
        if (numClusters == 1) { // centered about origin
            for (d = 0; d < UL_DIMENSION; d++) {
                center_point[d] = 0;
            }
        } else {
            for (d = 0; d < UL_DIMENSION; d++) {
                center_point[d] = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * 2 - 1;
            }
        }

        // iterate through each dimension, perturbing each by a separate gaussian
        for (d = 0; d < UL_DIMENSION; d++) { 

            // initialize the new gaussian
            std::normal_distribution<double> distribution(0.0, variance);

            // generate new gaussian sample along this dimension
            unsigned long long int oldIndex, newIndex;
            for (i = 0; i < UL_CLUSTERSIZE; i++) {

                // account for the shuffling -> new index for the point
                oldIndex = (c * UL_CLUSTERSIZE) + i;
                newIndex = (unsigned long long int) indices[oldIndex];

                // sample the gaussian to add some perturbation
                float perturbation = distribution(generator);
                dataPointer[newIndex*UL_DIMENSION + d] = center_point[d] + perturbation;
            }
        }
    }

    return;
}



// generate dataset of [N,D] floating point numbers
// uniformly distributed from [-1,1]
void sphericalDataset(float*& dataPointer, unsigned int dimension, unsigned int datasetSize, unsigned int seed) {
    float variance = 1; // variance of the gaussians

    // random number generator
    srand(seed);  // for same initializations of random variables
    std::default_random_engine generator;

    unsigned long long int numElements = (unsigned long long int)dimension * datasetSize;
    dataPointer = new float[numElements];

    // generate the values along each dimension with a new gaussian
    unsigned long long int UL_DIMENSION = (unsigned long long int) dimension;
    unsigned long long int d,i;
    for (d = 0; d < UL_DIMENSION; d++) {

        // new gaussian for each dimension
        std::normal_distribution<float> distribution(0.0, variance);

        // each point is a new sample of the gaussian
        for (i = 0; i < (unsigned long long int) datasetSize; i++) {
            dataPointer[i*UL_DIMENSION + d] = distribution(generator);
        }
    }

    //> Normalize Each Vector to Lie on the Sphere
    for (i = 0; i < (unsigned long long int) datasetSize; i++) {

        // calculate the magntiude of each vector
        float magnitude = 0;
        for (d = 0; d < dimension; d++) {
            magnitude += dataPointer[i*UL_DIMENSION + d]*dataPointer[i*UL_DIMENSION + d];
        }
        magnitude = sqrtf(magnitude);

        // divide each value by the magnitude
        for (d = 0; d < UL_DIMENSION; d++) {
            dataPointer[i*UL_DIMENSION + d] = dataPointer[i*UL_DIMENSION + d] / magnitude;
        }
    }

    return;
}

/**
 * @brief Datasets may not be randomized, i.e. ordered by class. Need to shuffle for optimization
 *
 * @param dataPointer
 * @param dimension
 * @param datasetSize
 */
void datasetShuffle(float*& dataPointer, unsigned int const dimension, unsigned int const datasetSize) {
    int seed = 7;
    std::mt19937 rng(seed);

    // create a vector for random assignment
    std::vector<unsigned int> indices(datasetSize, 0);
    std::iota(indices.begin(), indices.end(), 0);       // 0 , 1, 2, ... datasetSize-1
    std::shuffle(indices.begin(), indices.end(), rng);  // shuffle indices by rng seed

    unsigned long long int numElements = (unsigned long long int) datasetSize * (unsigned long long int) dimension;

    // copy dataPointer to different one for temp storage
    float* oldDataPointer = new float[numElements];
    for (unsigned long long int i = 0; i < numElements; i++) {
        oldDataPointer[i] = dataPointer[i];
    }

    // for all points in newly shuffled order
    unsigned long long int UL_DIMENSION = (unsigned long long int) dimension;
    unsigned long long int oldIndex, newIndex;
    unsigned long long int oldElement, newElement;
    for (newIndex = 0; newIndex < (unsigned long long int) datasetSize; newIndex++) {
        oldIndex = (unsigned long long int) indices[ (unsigned int) newIndex];

        // get all elements of the points
        for (unsigned long long int d = 0; d < UL_DIMENSION; d++) {
            newElement = (newIndex * UL_DIMENSION) + d;
            oldElement = (oldIndex * UL_DIMENSION) + d;
            dataPointer[newElement] = oldDataPointer[oldElement];
        }
    }

    // delete old
    delete[] oldDataPointer;
    return;
}

/**
 * @brief Load the Dataset and Testset onto one Pointer
 * 
 * @param data_path         // location of LA/ folder containing LA.txt and LA_query.txt
 * @param dataPointer 
 * @param dimension 
 * @param datasetSize 
 */
void LA(std::string data_path, float*& dataPointer, unsigned int& dimension, unsigned int& totalSize) {
    std::string datasetPath = std::string(data_path).append("LA.txt");
    std::string testsetPath = std::string(data_path).append("LA_query.txt");

    // hard-coded dataset parameters 
    dimension = 2;
    totalSize = 1073727 + 100; // N=1073727 in dataset, 100 in testset

    // putting the full set onto the float pointer
    dataPointer = new float[totalSize * dimension];
    unsigned int index = 0; // counter for the datapoints

    // open and save the data set
    // http://dbl.zju.edu.cn/~yjgao/MetricIndexes/Data.html
    std::ifstream file1(datasetPath.c_str());
    if (file1.is_open()) {
        std::string line,word;

        // first row as dimension, dataset size, distance metric
        std::getline(file1,line); 
        std::stringstream line_SS(line);
        getline(line_SS,word,' ');
        unsigned int dimension1 = (unsigned int) (std::atoi(word.c_str()));
        getline(line_SS,word,' ');
        unsigned int datasetSize1 = (unsigned int) (std::atoi(word.c_str()));
        getline(line_SS,word,' '); // euclidean distance
        printf("  * reading dataset file: N=%u, D=%u\n",datasetSize1, dimension1);

        // now, get all data
        while (std::getline(file1, line)) {
            if (line == "") break;

            // get the first dim
            std::stringstream lineSS1(line);
            std::getline(lineSS1, word, ' ');
            dataPointer[index * dimension + (0)] = std::atof(word.c_str());

            // get the second dim
            std::getline(lineSS1, word, ' ');
            dataPointer[index * dimension + (1)] = std::atof(word.c_str());
            index++;
        }
    }
    file1.close();

    // open and save the test set
    std::ifstream file2(testsetPath.c_str());
    if (file2.is_open()) {
        std::string line,word;

        // now, get all data
        while (std::getline(file2, line)) {
            if (line == "") break;

            // get the first dim
            std::stringstream lineSS1(line);
            std::getline(lineSS1, word, ' ');
            dataPointer[index * dimension + (0)] = std::atof(word.c_str());

            // get the second dim
            std::getline(lineSS1, word, ' ');
            dataPointer[index * dimension + (1)] = std::atof(word.c_str());
            index++;
        }
    }
    file2.close();
    printf("LA Dataset Opened. D=%u, N=%u\n",dimension, totalSize);

    return;
}

};  // namespace Datasets


#endif  // datasets_hpp