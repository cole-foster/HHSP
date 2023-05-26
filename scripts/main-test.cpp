// Cole Foster
// 2022-11-09
// Generalization of the HSP Graph
#include <chrono>
#include <cstdio>
#include <vector>

#include "CLI11.hpp"
#include "GHSP.hpp"
#include "NNS.hpp"
#include "datasets.hpp"
#include "pivot-index.hpp"

void printSet(std::vector<unsigned int> const &set);

int main(int argc, char **argv) {
    // printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataDirectory = "/users/cfoste18/scratch/sisap_data/";
    std::string dataset = "uniform";
    unsigned int dimension = 2;
    unsigned int datasetSize = 1000;
    unsigned int testsetSize = 100;
    std::vector<float> radiusVector{};
    int numThreads = 1;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size");
    app.add_option("-T,--testsetSize", testsetSize, "Testset Size");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-r,--radiusVector", radiusVector, "radius for the constant radius pivot index");
    CLI11_PARSE(app, argc, argv);


    // retrieve the dataset
    unsigned int totalSetSize = datasetSize + testsetSize;
    float* dataPointer;
    std::string metric = "euclidean";
    if (dataset == "uniform") {
        metric = "euclidean";
        Datasets::uniformDataset(dataPointer, dimension, totalSetSize, 3);
        datasetSize = totalSetSize - testsetSize; // incase it changes
    }  else if (dataset == "LA") { 
        std::string data_path = "/users/cfoste18/scratch/LA/";
        Datasets::LA(data_path, dataPointer, dimension, totalSetSize);
        if (datasetSize + testsetSize > totalSetSize) {
            datasetSize = totalSetSize - testsetSize;
        }
    } else {
        printf("Uknown dataset: %s. Aborting\n", dataset.c_str());
        return 0;
    }
    printf("d=%s,N=%u,D=%u\n", dataset.c_str(), datasetSize, dimension);
    unsigned long long int dStart, dEnd;
    std::chrono::high_resolution_clock::time_point tStart, tEnd;

    // create sparsematrix, really just holds datapointer and computes distances
    std::shared_ptr<SparseMatrix> sparseMatrix =
        std::make_shared<SparseMatrix>(dataPointer, totalSetSize, dimension);
    sparseMatrix->_datasetSize = datasetSize; 

    // create pivot index
    std::vector<PivotLayer> coverTree;
    PivotIndex::CoverTree_Greedy(radiusVector, *sparseMatrix, coverTree);
    printf("Cover Tree Valid: %u\n", PivotIndex::validateCoverTree(coverTree, *sparseMatrix));

    // perform hsp search by pivot index
    printf("HSP Search By PivotIndex: \n");
    std::vector<std::vector<unsigned int>> results_hsp_pivot{};
    results_hsp_pivot.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::GHSP_3L(datasetSize + queryIndex, coverTree, *sparseMatrix, results_hsp_pivot[queryIndex]);
    }
    tEnd = std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_hsp_pivot = (dEnd - dStart) / ((double)testsetSize);
    double time_hsp_pivot =
        std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double)testsetSize);
    printf("    * Time (ms): %.4f \n", time_hsp_pivot * 1000);
    printf("    * Distances: %.2f \n", distances_hsp_pivot);

    // perform hsp search by brute force
    printf("HSP Search By Brute Force: \n");
    std::vector<std::vector<unsigned int>> results_hsp_bf{};
    results_hsp_bf.resize(testsetSize);
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::HSP(datasetSize + queryIndex, datasetSize, *sparseMatrix, results_hsp_bf[queryIndex]);
    }

    // ensure correctness
    bool hsp_correct = true;
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        if (results_hsp_bf[queryIndex] != results_hsp_pivot[queryIndex]) {
            hsp_correct = false;
            printf("Error: Pivot NNS does not match GT!\n");
            printf("Q:%u,\n",queryIndex);
            printf("GT: "); printSet(results_hsp_bf[queryIndex]);
            printf("PI: "); printSet(results_hsp_pivot[queryIndex]);
            break;
        }
    }
    printf("    * Correct?: %u\n\n",hsp_correct);

    printf("Done! Have a good day! \n");
    delete[] dataPointer;
    return 0;
}

void printSet(std::vector<unsigned int> const &set) {
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u,", (*it1));
    }
    printf("}\n");
}