// Cole Foster
// November 29th, 2021
#include "datasets.hpp"

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
#include <sstream>
#include <fstream>

// generate dataset of [N,D] floating point numbers
// uniformly distributed from [-1,1]
void Datasets::uniformDataset(float*& dataPointer, unsigned int dimension, unsigned int N, float b, unsigned int seed) {
    srand(seed);  // for same initializations of random variables

    unsigned int numElements = (unsigned int)dimension * N;
    // delete[] dataPointer;
    dataPointer = new float[N * dimension];

    // produces numbers 0 to 1. multiply by 2b, subtract by b.
    for (int i = 0; i < numElements; i++) {
        dataPointer[i] = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX))*2*b - b;
    }
}


void Datasets::clusterDataset(float*& dataPointer, unsigned int dimension, unsigned int& datasetSize, std::string data_directory) {
    std::string filename = std::string(data_directory).append("cluster_").append(std::to_string(dimension)).append("D_N-1638400_single.bin");
    printf("Cluster Filename: %s\n",filename.c_str());

    // open dataset file
    std::ifstream inputFileStream(filename.c_str(), std::ios::binary);
    if(!inputFileStream.is_open()) {
        std::cout << "Error: cannot open file " << filename.c_str() << std::endl;
        return;
    }

    // read first integer, which should be the same as the dimension
    unsigned int data_dimension;
    inputFileStream.read((char*)&data_dimension, sizeof(int)); // reads first int, which gives us dimension of the following point
    inputFileStream.seekg(0, std::ios::end);

    // return error if dimension of dataset is different from input
    if (dimension != data_dimension) {
        std::cout << "Error: Dataset dimension different from input dimension" << std::endl;
        std::cout << "  - Input dimension: " << dimension << std::endl;
        std::cout << "  - Data dimension: " << data_dimension << std::endl;
        return;
    }

    // get number of points in the dataset
    std::ios::pos_type streamPosition = inputFileStream.tellg();
    std::size_t byteSize = (std::size_t)streamPosition;
    unsigned int fileDatasetSize = (unsigned int)(byteSize / (dimension*sizeof(float) + sizeof(int)));
    datasetSize = fileDatasetSize;

    // // return error if dataset size requested is larger than what we have
    // if (datasetSize > fileDatasetSize) {
    //     std::cout << "Error: Requested dataset size larger than dataset contains" << std::endl;
    //     std::cout << "  - Requested size: " << datasetSize << std::endl;
    //     std::cout << "  - True size: " << fileDatasetSize << std::endl;
    //     datasetSize = fileDatasetSize;
    // }
    unsigned int numElements = (unsigned int)dimension * datasetSize;
    dataPointer = new float[numElements];

    // read the elements of the dataset
    inputFileStream.seekg(0, std::ios::beg);
    for(std::size_t i = 0; i < datasetSize; i++)
    {
        inputFileStream.seekg(sizeof(int), std::ios::cur);
        inputFileStream.read((char*)(dataPointer + i*dimension), dimension*sizeof(float));
    }
    inputFileStream.close();
    return;
}


/**
 * @brief Datasets may not be randomized, i.e. ordered by class. Need to shuffle for optimization
 *
 * @param dataPointer
 * @param dimension
 * @param datasetSize
 */
void Datasets::datasetShuffle(float*& dataPointer, unsigned int dimension, unsigned int datasetSize) {
    int seed = 7;
    std::mt19937 rng(seed);

    std::vector<unsigned int> indices(datasetSize, 0);
    std::iota(indices.begin(), indices.end(), 0);       // 0 , 1, 2, ... datasetSize-1
    std::shuffle(indices.begin(), indices.end(), rng);  // shuffle indices by rng seed

    unsigned int numElements = datasetSize * dimension;

    // copy dataPointer to different one for temp storage
    float* oldDataPointer = new float[numElements];
    for (unsigned int i = 0; i < numElements; i++) {
        oldDataPointer[i] = dataPointer[i];
    }

    // for all points in newly shuffled order
    unsigned int oldIndex = 0, newIndex = 0, oldElement = 0, newElement = 0;
    for (newIndex = 0; newIndex < datasetSize; newIndex++) {
        unsigned int oldIndex = indices[newIndex];

        // get all elements of the points
        for (unsigned int d = 0; d < dimension; d++) {
            newElement = (newIndex * dimension) + d;
            oldElement = (oldIndex * dimension) + d;
            dataPointer[newElement] = oldDataPointer[oldElement];
        }
    }

    // delete old
    delete[] oldDataPointer;
    return;
}

/**
 * @brief Extract last few points from dataset as testset
 *
 * @param dimension
 * @param dataPointer
 * @param datasetSize
 * @param testPointer
 * @param testsetSize
 */
void Datasets::extractTestset(unsigned int dimension, float*& dataPointer, unsigned int& datasetSize, float*& testPointer, unsigned int testsetSize) {
    unsigned int numElements = (unsigned int)dimension * testsetSize;
    testPointer = new float[testsetSize * dimension];

    datasetSize -= testsetSize;
    unsigned int dataIndex = 0, testIndex = 0;
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        for (unsigned int d = 0; d < dimension; d++) {
            testIndex = (queryIndex)*dimension + d;
            dataIndex = (datasetSize + queryIndex) * dimension + d;
            testPointer[testIndex] = dataPointer[dataIndex];
        }
    }

    return;
}

/**
 * @brief Retrieve the WDBC dataset. 569 instances in 30D.
 *
 * @param dataPointer
 * @param dimension
 * @param datasetSize
 * @param dataDirectory
 */
void Datasets::Cities(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
    std::string data_path = std::string(dataDirectory).append("us_cities.csv");
    dimension = 2;
    datasetSize = 29880;
    dataPointer = new float[datasetSize * dimension];

    std::string line, word;
    std::fstream file(data_path.c_str(), std::ios::in);
    if (!file.is_open()) {
        printf("Open Error! \n");
    }
    if (!file.good()) {
        printf("not good!");
        return;
    }

    // get each line.
    std::getline(file, line); // heading

    unsigned int index = 0;
    while (std::getline(file, line)) {
        if (line == "") break;
        std::stringstream lineSS(line);

        // ignore first  values in each line. id 
        std::getline(lineSS, word, ',');
        for (int d = 0; d < dimension; d++) {
            std::getline(lineSS, word, ',');
            dataPointer[index * dimension + d] = std::atof(word.c_str());
        }
        index += 1;
    }

    return;
}

void Datasets::LA(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, 
           float*& testPointer, unsigned int& testSetSize, 
           std::string dataDirectory) {
    std::string trainPath = std::string(dataDirectory).append("LA.txt");
    std::string testPath = std::string(dataDirectory).append("LA_query.txt");

    // http://dbl.zju.edu.cn/~yjgao/MetricIndexes/Data.html
    std::ifstream file1(trainPath.c_str());
    if (file1.is_open()) {
        std::string line,word;

        // first row as dimension, dataset size, distance metric
        std::getline(file1,line); 
        std::stringstream line_SS(line);
        getline(line_SS,word,' ');
        dimension = (unsigned int) (std::atoi(word.c_str()));
        getline(line_SS,word,' ');
        datasetSize = (unsigned int) (std::atoi(word.c_str()));
        getline(line_SS,word,' '); // euclidean distance

        // now, get all data
        unsigned int index = 0;
        dataPointer = new float[datasetSize * dimension];
        while (std::getline(file1, line)) {
            if (line == "") break;

            // get the first dim
            std::stringstream lineSS1(line);
            std::getline(lineSS1, word, ' ');
            dataPointer[index * dimension + (0)] = std::atof(word.c_str());

            // get the second dim
            std::getline(lineSS1, word, ' ');
            dataPointer[index * dimension + (1)] = std::atof(word.c_str());
            index += 1;
        }
    }
    file1.close();

// 6032.63 5585.62
// 6033.83 5585.80
// 6026.04 5580.84
// 6044.37 5573.30
// 6042.56 5569.66

    std::ifstream file2(testPath.c_str());
    if (file2.is_open()) {
        std::string line,word;
        testSetSize = 100;

        // now, get all data
        unsigned int index = 0;
        testPointer = new float[testSetSize * dimension];
        while (std::getline(file2, line)) {
            if (line == "") break;

            // get the first dim
            std::stringstream lineSS1(line);
            std::getline(lineSS1, word, ' ');
            testPointer[index * dimension + (0)] = std::atof(word.c_str());

            // get the second dim
            std::getline(lineSS1, word, ' ');
            testPointer[index * dimension + (1)] = std::atof(word.c_str());
            index += 1;
        }
    }
    file2.close();
    return;
}





// void SIFT1M(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string data_path) {

//     printf("get SIFT1M...\n");

//     // from https://github.com/facebookresearch/faiss/blob/main/demos/demo_sift1M.cpp
//     FILE* f = fopen(data_path.c_str(), "r");
//     if (!f) {
//         fprintf(stderr, "could not open %s\n", data_path.c_str());
//         perror("");
//         abort();
//     }

//     printf("opened...\n");
//     int d;
//     fread(&d, 1, sizeof(int), f);
//     //assert((d > 0 && d < 1000000) || !"unreasonable dimension");
//     fseek(f, 0, SEEK_SET);
//     struct stat st;
//     fstat(fileno(f), &st);
//     size_t sz = st.st_size;
//     //assert(sz % ((d + 1) * 4) == 0 || !"weird file size");
//     size_t n = sz / ((d + 1) * 4);

//     dataPointer = new float[n * (d + 1)];
//     size_t nr = fread(dataPointer, sizeof(float), n * (d + 1), f);
//     //assert(nr == n * (d + 1) || !"could not read whole file");

//     // shift array to remove row headers
//     for (size_t i = 0; i < n; i++) {
//         memmove(dataPointer + i * d, dataPointer + 1 + i * (d + 1), d * sizeof(*dataPointer));
//     }

//     printf("done!\n");

//     dimension = (unsigned int) d;
//     datasetSize = (unsigned int) n;
//     fclose(f);
//     return;
// }


// void DEEP10M(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
//     std::string data_path = std::string(dataDirectory).append("deep10M.fvecs");
//     //deep10M.fvecs

//     printf("get DEEP...\n");

//     // from https://github.com/facebookresearch/faiss/blob/main/demos/demo_sift1M.cpp
//     FILE* f = fopen(data_path.c_str(), "r");
//     if (!f) {
//         fprintf(stderr, "could not open %s\n", data_path.c_str());
//         perror("");
//         abort();
//     }

//     printf("opened...\n");
//     int d;
//     fread(&d, 1, sizeof(int), f);
//     //assert((d > 0 && d < 1000000) || !"unreasonable dimension");
//     fseek(f, 0, SEEK_SET);
//     struct stat st;
//     fstat(fileno(f), &st);
//     size_t sz = st.st_size;
//     //assert(sz % ((d + 1) * 4) == 0 || !"weird file size");
//     size_t n = sz / ((d + 1) * 4);

//     //if (datasetSize > 0) {
//     //    n = (size_t) (datasetSize);
//     //}

//     dataPointer = new float[n * (d + 1)];
//     size_t nr = fread(dataPointer, sizeof(float), n * (d + 1), f);
//     //assert(nr == n * (d + 1) || !"could not read whole file");

//     // shift array to remove row headers
//     for (size_t i = 0; i < n; i++) {
//         memmove(dataPointer + i * d, dataPointer + 1 + i * (d + 1), d * sizeof(*dataPointer));
//     }

//     printf("done!\n");

//     dimension = (unsigned int) d;
//     datasetSize = (unsigned int) n;
//     fclose(f);
//     return;
// }