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
#include "animation.hpp"

void printSet(std::vector<unsigned int> const &set);

int main(int argc, char **argv) {
    // printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataset = "uniform";
    unsigned int dimension = 2;
    float cube_length = 1;
    unsigned int datasetSize = 1000;
    unsigned int testsetSize = 100;
    std::vector<float> radiusVector{};
    int numThreads = 1;
    bool verbose = true;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-c,--cube_length", cube_length, "Length of Cube Side");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size");
    app.add_option("-T,--testsetSize", testsetSize, "Testset Size");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-r,--radius", radiusVector, "radius for the constant radius pivot index");
    app.add_option("-v,--verbose", verbose, "Show hRNG Statistics During Eval");
    CLI11_PARSE(app, argc, argv);

    //====================================================================
    //                      Getting Dataset
    //====================================================================

    // get dataset (with test set added)
    // dataset: 0...datsetSize-1, queryset: datasetSize-datasetSize+testsetSize
    float *dataPointer = NULL;
    if (dataset == "uniform") {
        Datasets::uniformDataset(dataPointer, dimension, datasetSize + testsetSize, cube_length, 3);
    } else {
        printf("Unrecognized dataset: %s\n", dataset.c_str());
        return 0;
    }

    unsigned long long int dStart, dEnd;
    std::chrono::high_resolution_clock::time_point tStart, tEnd;

    //====================================================================
    //                Initializing Distance Data Structure
    //====================================================================

    // create sparsematrix, really just holds datapointer and computes distances
    std::shared_ptr<SparseMatrix> sparseMatrix =
        std::make_shared<SparseMatrix>(dataPointer, datasetSize + testsetSize, dimension);
    sparseMatrix->_datasetSize = datasetSize;

    //====================================================================
    //                     Creating Pivot Index
    //====================================================================
    printf("Constructing the Pivot-Index: \n");
    int numberOfLayers = radiusVector.size() + 1;
    std::vector<Pivot> pivotsList;
    std::vector<unsigned int> pivotsPerLayer{};
    PivotIndex::Greedy_MultiLayer(radiusVector, *sparseMatrix, pivotsList, pivotsPerLayer);
    pivotsPerLayer.push_back(datasetSize);
    printf("    * PivotsPerLayer: \n");
    printSet(pivotsPerLayer);

    //====================================================================
    //                      HSP Search
    //====================================================================

    // perform hsp search by pivot index
    printf("HSP Search By PivotIndex: \n");
    std::vector<std::vector<unsigned int>> results_hsp_pivot{};
    results_hsp_pivot.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::GHSP_Search(datasetSize + queryIndex, pivotsList, *sparseMatrix, results_hsp_pivot[queryIndex]);
    }
    tEnd = std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_hsp_pivot = (dEnd - dStart) / ((double)testsetSize);
    double time_hsp_pivot =
        std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double)testsetSize);
    printf("    * Time (ms): %.4f \n", time_hsp_pivot * 1000);
    printf("    * Distances: %.2f \n", distances_hsp_pivot);


    //====================================================================
    //                      Animation
    //====================================================================



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
