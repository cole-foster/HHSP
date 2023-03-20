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
// generate dataset of [N,D] floating point numbers
// uniformly distributed from [-1,1]
void uniformDataset(float*& dataPointer, unsigned int dimension, unsigned int N, unsigned int seed) {
    srand(seed);  // for same initializations of random variables

    unsigned int numElements = (unsigned int)dimension * N;
    // delete[] dataPointer;
    dataPointer = new float[N * dimension];

    for (int i = 0; i < numElements; i++) {
        dataPointer[i] = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * 2 - 1;
    }
}

void clusterDataset(float*& dataPointer, unsigned int dimension, unsigned int& datasetSize, std::string data_directory) {
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
void datasetShuffle(float*& dataPointer, unsigned int dimension, unsigned int datasetSize) {
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
void extractTestset(unsigned int dimension, float*& dataPointer, unsigned int& datasetSize, float*& testPointer, unsigned int testsetSize) {
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
void Cities(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
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

/**
 * @brief Retrieve the Corel68k dataset. 68040 instances in 57D. Takes Color Histograms, Color Moments, and Color Texture
 *
 * @param dataPointer
 * @param dimension
 * @param datasetSize
 * @param dataDirectory
 */
void Corel68k1(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
    std::string data_path_H = std::string(dataDirectory).append("ColorHistogram.asc");
    dimension = 32;
    datasetSize = 68040;
    dataPointer = new float[datasetSize * dimension];

    std::string line_H, word_H;
    std::fstream file_H(data_path_H.c_str(), std::ios::in);
    if (!file_H.is_open()) {
        printf("Open Error! \n");
    }

    // get each line.
    unsigned int index = 0;
    while (std::getline(file_H, line_H)) {
        if (line_H == "") break;

        // get the 32D color histogram first.
        unsigned int sd = 0;
        std::stringstream lineSS_H(line_H);
        std::getline(lineSS_H, word_H, ' ');  // skip first entry as id
        for (int d = 0; d < 32; d++) {
            std::getline(lineSS_H, word_H, ' ');
            dataPointer[index*dimension + (sd + d)] = std::atof(word_H.c_str());
        }

        index += 1;
    }

    return;
}

/**
 * @brief Retrieve the Corel68k dataset. 68040 instances in 57D. Takes Color Histograms, Color Moments, and Color Texture
 *
 * @param dataPointer
 * @param dimension
 * @param datasetSize
 * @param dataDirectory
 */
void Corel68k2(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
    std::string data_path_M = std::string(dataDirectory).append("ColorMoments.asc");
    dimension = 9;
    datasetSize = 68040;
    dataPointer = new float[datasetSize * dimension];

    std::string line_M, word_M;
    std::fstream file_M(data_path_M.c_str(), std::ios::in);
    if (!file_M.is_open()) {
        printf("Open Error! \n");
    }

    // get each line.
    unsigned int index = 0;
    while (std::getline(file_M, line_M)) {
        if (line_M == "") break;

        // get the 9D color moment next.
        unsigned int sd = 0;
        std::stringstream lineSS_M(line_M);
        std::getline(lineSS_M, word_M, ' ');  // skip first entry as id
        for (int d = 0; d < 9; d++) {
            std::getline(lineSS_M, word_M, ' ');
            dataPointer[index * dimension + (sd + d)] = std::atof(word_M.c_str());
        }

        index += 1;
    }

    return;
}

/**
 * @brief Retrieve the Corel68k dataset. 68040 instances in 57D. Takes Color Histograms, Color Moments, and Color Texture
 *
 * @param dataPointer
 * @param dimension
 * @param datasetSize
 * @param dataDirectory
 */
void Corel68k3(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
    std::string data_path_T = std::string(dataDirectory).append("CoocTexture.asc");
    dimension = 16;
    datasetSize = 68040;
    dataPointer = new float[datasetSize * dimension];

    std::string line_T, word_T;
    std::fstream file_T(data_path_T.c_str(), std::ios::in);
    if (!file_T.is_open()) {
        printf("Open Error! \n");
    }

    // get each line.
    unsigned int index = 0;
    while (std::getline(file_T, line_T)) {
        if (line_T == "") break;

        // get the 16D co-occurence texture next.
        unsigned int sd = 0;
        std::stringstream lineSS_T(line_T);
        std::getline(lineSS_T, word_T, ' ');  //  skip first entry as id
        for (int d = 0; d < 16; d++) {
            std::getline(lineSS_T, word_T, ' ');
            dataPointer[index * dimension + (sd + d)] = std::atof(word_T.c_str());
        }

        index += 1;
    }

    return;
}


void SIFT1M(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string data_path) {

    printf("get SIFT1M...\n");

    // from https://github.com/facebookresearch/faiss/blob/main/demos/demo_sift1M.cpp
    FILE* f = fopen(data_path.c_str(), "r");
    if (!f) {
        fprintf(stderr, "could not open %s\n", data_path.c_str());
        perror("");
        abort();
    }

    printf("opened...\n");
    int d;
    fread(&d, 1, sizeof(int), f);
    //assert((d > 0 && d < 1000000) || !"unreasonable dimension");
    fseek(f, 0, SEEK_SET);
    struct stat st;
    fstat(fileno(f), &st);
    size_t sz = st.st_size;
    //assert(sz % ((d + 1) * 4) == 0 || !"weird file size");
    size_t n = sz / ((d + 1) * 4);

    dataPointer = new float[n * (d + 1)];
    size_t nr = fread(dataPointer, sizeof(float), n * (d + 1), f);
    //assert(nr == n * (d + 1) || !"could not read whole file");

    // shift array to remove row headers
    for (size_t i = 0; i < n; i++) {
        memmove(dataPointer + i * d, dataPointer + 1 + i * (d + 1), d * sizeof(*dataPointer));
    }

    printf("done!\n");

    dimension = (unsigned int) d;
    datasetSize = (unsigned int) n;
    fclose(f);
    return;
}


void DEEP10M(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
    std::string data_path = std::string(dataDirectory).append("deep10M.fvecs");
    //deep10M.fvecs

    printf("get DEEP...\n");

    // from https://github.com/facebookresearch/faiss/blob/main/demos/demo_sift1M.cpp
    FILE* f = fopen(data_path.c_str(), "r");
    if (!f) {
        fprintf(stderr, "could not open %s\n", data_path.c_str());
        perror("");
        abort();
    }

    printf("opened...\n");
    int d;
    fread(&d, 1, sizeof(int), f);
    //assert((d > 0 && d < 1000000) || !"unreasonable dimension");
    fseek(f, 0, SEEK_SET);
    struct stat st;
    fstat(fileno(f), &st);
    size_t sz = st.st_size;
    //assert(sz % ((d + 1) * 4) == 0 || !"weird file size");
    size_t n = sz / ((d + 1) * 4);

    //if (datasetSize > 0) {
    //    n = (size_t) (datasetSize);
    //}

    dataPointer = new float[n * (d + 1)];
    size_t nr = fread(dataPointer, sizeof(float), n * (d + 1), f);
    //assert(nr == n * (d + 1) || !"could not read whole file");

    // shift array to remove row headers
    for (size_t i = 0; i < n; i++) {
        memmove(dataPointer + i * d, dataPointer + 1 + i * (d + 1), d * sizeof(*dataPointer));
    }

    printf("done!\n");

    dimension = (unsigned int) d;
    datasetSize = (unsigned int) n;
    fclose(f);
    return;
}


// https://stackoverflow.com/questions/8286668/how-to-read-mnist-data-in-c
int reverseInt (int i) 
{
    unsigned char c1, c2, c3, c4;

    c1 = i & 255;
    c2 = (i >> 8) & 255;
    c3 = (i >> 16) & 255;
    c4 = (i >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
}

void MNIST(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, std::string dataDirectory) {
    std::string trainPath = std::string(dataDirectory).append("mnist-train.ubyte");
    std::string testPath = std::string(dataDirectory).append("mnist-test.ubyte");

    // https://deepai.org/dataset/mnist
    std::ifstream file1(trainPath.c_str());
    unsigned int datasetSize1 = 0;
    if (file1.is_open()) {
        int magic_number=0;
        int number_of_images=0;
        int n_rows=0;
        int n_cols=0;
        file1.read((char*)&magic_number,sizeof(magic_number)); 
        magic_number= reverseInt(magic_number);
        file1.read((char*)&number_of_images,sizeof(number_of_images));
        number_of_images= reverseInt(number_of_images);
        file1.read((char*)&n_rows,sizeof(n_rows));
        n_rows= reverseInt(n_rows);
        file1.read((char*)&n_cols,sizeof(n_cols));
        n_cols= reverseInt(n_cols);
        dimension = (unsigned int)(n_rows*n_cols);
        datasetSize1 = (unsigned int) number_of_images;
    }

    // https://deepai.org/dataset/mnist
    unsigned int datasetSize2 = 0;
    std::ifstream file2(testPath.c_str());
    if (file2.is_open()) {
        int magic_number=0;
        int number_of_images=0;
        int n_rows=0;
        int n_cols=0;
        file2.read((char*)&magic_number,sizeof(magic_number)); 
        magic_number= reverseInt(magic_number);
        file2.read((char*)&number_of_images,sizeof(number_of_images));
        number_of_images= reverseInt(number_of_images);
        file2.read((char*)&n_rows,sizeof(n_rows));
        n_rows= reverseInt(n_rows);
        file2.read((char*)&n_cols,sizeof(n_cols));
        n_cols= reverseInt(n_cols);
        datasetSize2 = (unsigned int) number_of_images;
    }

    // now read all
    datasetSize = datasetSize1 + datasetSize2;
    dataPointer = new float[datasetSize * (dimension)];
    if (file1.is_open()) {
        for (unsigned int i = 0; i < datasetSize1; ++i) {
            for (unsigned int d = 0; d < dimension; ++d) {
                unsigned char temp = 0;
                file1.read((char*)&temp,sizeof(temp));
                dataPointer[(i*dimension)+d] = (float)(temp);
            }
        }
    }
    unsigned int start_idx = datasetSize1*dimension; 
    if (file2.is_open()) {
        for (unsigned int i = 0; i < datasetSize2; ++i) {
            for (unsigned int d = 0; d < dimension; ++d) {
                unsigned char temp = 0;
                file2.read((char*)&temp,sizeof(temp));
                dataPointer[start_idx+(i*dimension)+d] = (float)(temp);
            }
        }
    }

    file1.close();
    file2.close();
    return;
}

void LA(float*& dataPointer, unsigned int& dimension, unsigned int& datasetSize, 
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


};  // namespace Datasets



#endif  // datasets_hpp