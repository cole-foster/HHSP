// Cole Foster
// 2022-11-09
// Generalization of the HSP Graph
#include <chrono>
#include <cstdio>
#include <vector>

#include "CLI11.hpp"
#include "datasets.hpp"
#include "pivot-index.hpp"
#include "NNS.hpp"
#include "GHSP.hpp"

void printSet(std::vector<unsigned int> const &set);


int main(int argc, char **argv) {
    // printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataDirectory = "/users/cfoste18/scratch/sisap_data/";
    std::string dataset = "uniform";
    unsigned int dimension = 2;
    float cube_length = 1;
    unsigned int datasetSize = 1000;
    unsigned int testsetSize = 100;
    std::vector<float> radiusVector{};
    int numThreads = 1;
    bool verbose = true;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-c,--cube_length", cube_length, "cube_length");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size");
    app.add_option("-T,--testsetSize", testsetSize, "Testset Size");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-r,--radius", radiusVector, "radius for the constant radius pivot index");
    app.add_option("-v,--verbose", verbose, "Show hRNG Statistics During Eval");
    CLI11_PARSE(app, argc, argv);

    // get dataset (with test set added)
    // dataset: 0...datsetSize-1, queryset: datasetSize-datasetSize+testsetSize
    float *dataPointer = NULL;
    if (dataset == "uniform") {
        Datasets::uniformDataset(dataPointer, dimension, datasetSize + testsetSize,cube_length, 3);
    } else if (dataset == "sphere") {
        Datasets::sphericalDistribution(dataPointer, dimension, datasetSize + testsetSize,1,3);
    } else if (dataset == "SIFT") {
        unsigned int startSize = datasetSize;
        Datasets::SIFT1M(dataPointer, dimension, datasetSize);
        if (startSize < datasetSize - testsetSize) {
            datasetSize = startSize;   
        } else {
            printf("Not enough data. Requested N=%u points with T=%u queries, there is only %u\n",startSize, testsetSize, datasetSize);
            datasetSize = datasetSize - testsetSize;
        }
    } else {
        printf("Unrecognized dataset: %s\n", dataset.c_str());
        return 0;
    }

    unsigned long long int dStart, dEnd;
    std::chrono::high_resolution_clock::time_point tStart,tEnd;

    // create sparsematrix, really just holds datapointer and computes distances
    std::shared_ptr<SparseMatrix> sparseMatrix = std::make_shared<SparseMatrix>(dataPointer, datasetSize+testsetSize, dimension);
    sparseMatrix->_datasetSize = datasetSize;
    int numberOfLayers = radiusVector.size()+1;

    // create pivot index
    printf("Pivot Selection: \n");
    printf("    * N=%u\n",datasetSize);
    printf("    * D=%u\n",dimension);
    std::vector<PivotLayer> pivotLayers;    
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    if (numberOfLayers == 2) {
        PivotIndex::Greedy_2L(radiusVector[0], *sparseMatrix, pivotLayers);
    } else if (numberOfLayers == 3) {
        PivotIndex::Greedy_3L(radiusVector, *sparseMatrix, pivotLayers);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    unsigned long long int distances_ps = (dEnd - dStart);
    double time_ps = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);
    printf("    * |P|=%u\n",pivotLayers[0].pivotIndices->size());    
    printf("    * PS Distances: %llu \n",distances_ps);
    printf("    * PS Time (s): %.4f \n",time_ps);
    // bool success = false;
    // // if (numberOfLayers == 2) {
    // //     bool success = PivotSelection::validatePivotSelection(pivotLayer,sparseMatrix);
    // if (numberOfLayers == 3) {
    //     success = PivotIndex::validatePivotSelection_3L(pivotLayers, *sparseMatrix);
    // }
    // printf("    * Minimal Coverage?: %u\n\n",success);

    // perform nns time for all pivots, get average results
    printf("Nearest Neighbor Search By Brute Force: \n");
    std::vector<unsigned int> results_nns_bf{};
    results_nns_bf.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        NNS::NNS_BruteForce(datasetSize + queryIndex, datasetSize, *sparseMatrix, results_nns_bf[queryIndex]);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_nns_brute = (dEnd - dStart)/((double)testsetSize);
    double time_nns_brute = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);
    printf("    * Time (ms): %.4f \n",time_nns_brute*1000);
    printf("    * Distances: %.2f \n\n",distances_nns_brute);

    // perform nns by index
    printf("Nearest Neighbor Search By PivotIndex: \n");
    std::vector<unsigned int> results_nns_pivot{};
    results_nns_pivot.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    if (numberOfLayers == 2) {
        for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
            NNS::NNS_2L(datasetSize + queryIndex, pivotLayers, *sparseMatrix, results_nns_pivot[queryIndex]);
        }
    } else if (numberOfLayers == 3) {
        for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
            NNS::NNS_3L(datasetSize + queryIndex, pivotLayers, *sparseMatrix, results_nns_pivot[queryIndex]);
        }
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_nns_pivot= (dEnd - dStart)/((double)testsetSize);
    double time_nns_pivot = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);
    printf("    * Time (ms): %.4f \n",time_nns_pivot*1000);
    printf("    * Distances: %.2f \n",distances_nns_pivot);

    // ensure correctness
    bool nns_correct = true;
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        if (results_nns_bf[queryIndex] != results_nns_pivot[queryIndex]) {
            nns_correct = false;
            printf("Error: Pivot NNS does not match GT!\n");
            printf("Q:%u, GT:%u, Pivot: %u\n",queryIndex,results_nns_bf[queryIndex],results_nns_pivot[queryIndex]);
            break;
        }
    }
    printf("    * Correct?: %u\n\n",nns_correct);

    // perform hsp search by brute force
    printf("HSP Search By Brute Force: \n");
    std::vector<std::vector<unsigned int>> results_hsp_bf{};
    results_hsp_bf.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();\
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::HSP(datasetSize + queryIndex, datasetSize, *sparseMatrix, results_hsp_bf[queryIndex]);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_hsp_brute = (dEnd - dStart)/((double)testsetSize);
    double time_hsp_brute = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);
    printf("    * Time (ms): %.4f \n",time_hsp_brute*1000);
    printf("    * Distances: %.2f \n\n",distances_hsp_brute);

    // perform hsp search by pivot index
    printf("HSP Search By PivotIndex: \n");
    std::vector<std::vector<unsigned int>> results_hsp_pivot{};
    results_hsp_pivot.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    if (numberOfLayers == 2) {
        for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
            GHSP::GHSP_2L(datasetSize + queryIndex, pivotLayers[0], *sparseMatrix, results_hsp_pivot[queryIndex]);
        }
    } else if (numberOfLayers == 3) {
        for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
            GHSP::GHSP_3L(datasetSize + queryIndex, pivotLayers, *sparseMatrix, results_hsp_pivot[queryIndex]);
            // break;
        }
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_hsp_pivot = (dEnd - dStart)/((double)testsetSize);
    double time_hsp_pivot = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);
    printf("    * Time (ms): %.4f \n",time_hsp_pivot*1000);
    printf("    * Distances: %.2f \n",distances_hsp_pivot);

    // // ensure correctness
    double ave_hsp = 0;//(double) results_hsp_bf[0].size();
    bool hsp_correct = true;
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        ave_hsp += (double) results_hsp_bf[queryIndex].size();
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
    ave_hsp /= (double) testsetSize;

    printf("D,N,r,|P|,Index Distances,Index Time (s),BF NN Distances,BF NN Time (ms), Pivot NN Distances, Pivot NN Time (ms), BF HSP Distances, BF HSP Time (ms), Pivot HSP Distances, Pivot HSP Time (ms),Correct, Ave HSP\n");
    printf("%u,%u,%.6f,%u,",dimension,datasetSize,radiusVector[0],pivotLayers[0].pivotIndices->size());
    printf("%llu,%.4f,",distances_ps,time_ps);
    printf("%.2f,%.4f,",distances_nns_brute,time_nns_brute*1000);
    printf("%.2f,%.4f,",distances_nns_pivot,time_nns_pivot*1000);
    printf("%.2f,%.4f,",distances_hsp_brute,time_hsp_brute*1000);
    printf("-,-,");
    printf("%.2f,%.4f,",distances_hsp_pivot,time_hsp_pivot*1000);
    printf("%u,",hsp_correct);
    // printf("-,");
    printf("%.2f\n",ave_hsp);

    printf("Done! Have a good day! \n");
    delete[] dataPointer;
    return 0;
}


void printSet(std::vector<unsigned int> const &set)
{
    printf("{");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++)
    {
        printf("%u,", (*it1));
    }
    printf("}\n");
}


