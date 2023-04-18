// Cole Foster
// 2022-11-09
// Generalization of the HSP Graph
#include <chrono>
#include <cstdio>
#include <vector>

#include "CLI11.hpp"
#include "datasets.hpp"
#include "GHSP.hpp"
#include "pivot-selection.hpp"

void printSet(std::vector<unsigned int> const &set);

double functionEval(float const radius, SparseMatrix& sparseMatrix, unsigned int const testsetSize);
void fullRun(float const radius, SparseMatrix& sparseMatrix, unsigned int const testsetSize);

int main(int argc, char **argv) {
    // printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataDirectory = "/users/cfoste18/scratch/sisap_data/";
    std::string dataset = "uniform";
    unsigned int dimension = 2;
    float cube_length = 1;
    unsigned int datasetSize = 1000;
    unsigned int testsetSize = 100;
    int numThreads = 1;
    bool verbose = true;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-c,--cube_length", cube_length, "cube_length");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size");
    app.add_option("-T,--testsetSize", testsetSize, "Testset Size");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-v,--verbose", verbose, "Show hRNG Statistics During Eval");

    // optimization parameters
    float a = 0;
    float b = 1;
    float tol = 0.001;
    app.add_option("-a,--a", a, "left bound");
    app.add_option("-b,--b", b, "right bound");
    app.add_option("-t,--t", tol, "tolerance for stopping");
    CLI11_PARSE(app, argc, argv);

    // get dataset (with test set added)
    // dataset: 0...datsetSize-1, queryset: datasetSize-datasetSize+testsetSize
    float *dataPointer = NULL;
    if (dataset == "uniform") {
        Datasets::uniformDataset(dataPointer, dimension, datasetSize + testsetSize, cube_length, 3);
    } else if (dataset == "sphere") {
        Datasets::sphericalDistribution(dataPointer, dimension, datasetSize + testsetSize, 1, 3);
    } else if (dataset == "SIFT") {
        unsigned int startSize = datasetSize;
        Datasets::SIFT1M(dataPointer, dimension, datasetSize);
        Datasets::datasetShuffle(dataPointer, dimension, datasetSize);
        printf("SIFT Dataset: N=%u, D=%u\n",datasetSize, dimension);
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

    // create sparsematrix, really just holds datapointer and computes distances
    std::shared_ptr<SparseMatrix> sparseMatrix = std::make_shared<SparseMatrix>(dataPointer, datasetSize+testsetSize, dimension);
    sparseMatrix->_datasetSize = datasetSize;

    // optimization set up
    float R = (1+sqrt(5))/2; // golden ratio
    float sect_c = a + (b-a)/R; // left middle
    float sect_d = b - (b-a)/R; // right middle

    // one function eval per iteration
    bool flag_eval_c = true;
    double fsect_c;
    printf("D,N,r,|P|,Index Distances,Index Time (s), Pivot HSP Distances, Pivot HSP Time (ms)\n");
    double fsect_d = functionEval(sect_d,*sparseMatrix,testsetSize);

    // optimization loop
    while(fabs(a-b)/2 > tol) { 

        // one function evaluation
        if (flag_eval_c) {
            fsect_c = functionEval(sect_c,*sparseMatrix,testsetSize);
        } else {
            fsect_d = functionEval(sect_d,*sparseMatrix,testsetSize);
        }

        if (fsect_c < fsect_d) { // need to re-eval sect_c
            flag_eval_c = true;
            fsect_d = fsect_c; 

            a = sect_d;
            sect_d = sect_c;
            sect_c = a + (b-a)/R;
        } else { // need to re-eval sect_d
            flag_eval_c = false;
            fsect_c = fsect_d; 

            b = sect_c;
            sect_c = sect_d;
            sect_d = b - (b-a)/R;
        }
    }

    printf("--------------- Optimal ------------------\n");
    double optimalRadius = (a+b)/2;
    fullRun(optimalRadius,*sparseMatrix,testsetSize);

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


double functionEval(float const radius, SparseMatrix& sparseMatrix, unsigned int const testsetSize) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    unsigned int const dimension = sparseMatrix._dimension;
    unsigned long long int dStart, dEnd;
    std::chrono::high_resolution_clock::time_point tStart,tEnd;

    // create pivot index
    PivotLayer pivotLayer;    
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    PivotSelection::selectPivots(radius, sparseMatrix, pivotLayer);
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    unsigned long long int distances_ps = (dEnd - dStart);
    double time_ps = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();

    // perform nns by index
    unsigned int nns_neighbor;
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::pivotNNS(datasetSize + queryIndex, pivotLayer, sparseMatrix, nns_neighbor);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    double distances_nns_pivot= (dEnd - dStart)/((double)testsetSize);
    double time_nns_pivot = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);

    // perform hsp search by pivot index
    double ave_neighbors = 0;
    std::vector<unsigned int> neighbors{};
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::GHSP_Search(datasetSize + queryIndex, pivotLayer, sparseMatrix, neighbors);
        ave_neighbors += (double) neighbors.size();
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    double distances_hsp_pivot = (dEnd - dStart)/((double)testsetSize);
    double time_hsp_pivot = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);
    
    printf("%u,%u,%.6f,%u,",dimension,datasetSize,radius,pivotLayer.pivotIndices->size());
    printf("%llu,%.4f,",distances_ps,time_ps);
    printf("%.2f,%.4f,",distances_nns_pivot,time_nns_pivot*1000);
    printf("%.2f,%.4f,%.2f\n",distances_hsp_pivot,time_hsp_pivot*1000, ave_neighbors / ((double) testsetSize));
    return distances_nns_pivot; //distances_hsp_pivot;
}


void fullRun(float const radius, SparseMatrix& sparseMatrix, unsigned int const testsetSize) {
    unsigned int const datasetSize = sparseMatrix._datasetSize;
    unsigned int const dimension = sparseMatrix._dimension;
    unsigned long long int dStart, dEnd;
    std::chrono::high_resolution_clock::time_point tStart,tEnd;

    // create pivot index
    PivotLayer pivotLayer;    
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    PivotSelection::selectPivots(radius, sparseMatrix, pivotLayer);
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    unsigned long long int distances_ps = (dEnd - dStart);
    double time_ps = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();

    // perform nns time for all pivots, get average results
    std::vector<unsigned int> results_nns_bf{};
    results_nns_bf.resize(testsetSize);
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::bruteNNS(datasetSize + queryIndex, datasetSize, sparseMatrix, results_nns_bf[queryIndex]);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    double distances_nns_brute = (dEnd - dStart)/((double)testsetSize);
    double time_nns_brute = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);

    // perform nns by index
    std::vector<unsigned int> results_nns_pivot{};
    results_nns_pivot.resize(testsetSize);
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::pivotNNS(datasetSize + queryIndex, pivotLayer, sparseMatrix, results_nns_pivot[queryIndex]);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    double distances_nns_pivot= (dEnd - dStart)/((double)testsetSize);
    double time_nns_pivot = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);

    // perform hsp search by brute force
    std::vector<std::vector<unsigned int>> results_hsp_bf{};
    results_hsp_bf.resize(testsetSize);
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();\
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::HSP(datasetSize + queryIndex, datasetSize, sparseMatrix, results_hsp_bf[queryIndex]);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    double distances_hsp_brute = (dEnd - dStart)/((double)testsetSize);
    double time_hsp_brute = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);

    // perform hsp search by pivot index
    std::vector<std::vector<unsigned int>> results_hsp_pivot{};
    results_hsp_pivot.resize(testsetSize);
    dStart = sparseMatrix._distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();\
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::GHSP_Search(datasetSize + queryIndex, pivotLayer, sparseMatrix, results_hsp_pivot[queryIndex]);
    }
    tEnd= std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix._distanceComputationCount;
    double distances_hsp_pivot = (dEnd - dStart)/((double)testsetSize);
    double time_hsp_pivot = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double) testsetSize);

    // ensure correctness
    double ave_deg = 0;
    bool hsp_correct = true;
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        ave_deg += (double) results_hsp_bf[queryIndex].size();
        if (results_hsp_bf[queryIndex] != results_hsp_pivot[queryIndex]) {
            hsp_correct = false;
            break;
        }
    }

    printf("D,N,r,|P|,Index Distances,Index Time (s),BF NN Distances,BF NN Time (ms), Pivot NN Distances, Pivot NN Time (ms), BF HSP Distances, BF HSP Time (ms), Pivot HSP Distances, Pivot HSP Time (ms),Correct\n");
    printf("%u,%u,%.6f,%u,",dimension,datasetSize,radius,pivotLayer.pivotIndices->size());
    printf("%llu,%.4f,",distances_ps,time_ps);
    printf("%.2f,%.4f,",distances_nns_brute,time_nns_brute*1000);
    printf("%.2f,%.4f,",distances_nns_pivot,time_nns_pivot*1000);
    printf("%.2f,%.4f,",distances_hsp_brute,time_hsp_brute*1000);
    printf("%.2f,%.4f,",distances_hsp_pivot,time_hsp_pivot*1000);
    printf("%u,",hsp_correct);
    printf("%.2f\n",ave_deg/((double)testsetSize));

    return;
}