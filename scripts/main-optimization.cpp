// Cole Foster
// 2022-11-09
// Generalization of the HSP Graph
#include <chrono>
#include <cstdio>
#include <nlopt.hpp>
#include <vector>

#include "CLI11.hpp"
#include "datasets.hpp"
#include "GHSP.hpp"
#include "NNS.hpp"
#include "pivot-index.hpp"

struct GHSP_Data {
    GHSP_Data(){};
    GHSP_Data(unsigned int const dimension, float *const &dataPointer, unsigned int const datasetSize,
              unsigned int const testsetSize)
        : _dimension(dimension), _dataPointer(dataPointer), _datasetSize(datasetSize), _testsetSize(testsetSize){};
    ~GHSP_Data(){};

    // build and search parameters
    unsigned int _dimension = 0;
    float *_dataPointer = NULL;
    unsigned int _datasetSize = 0;
    unsigned int _testsetSize = 0;
};

double objectiveFunction(const std::vector<double> &x, std::vector<double> &grad, void *data_ghsp);
bool COBLYA(unsigned n, std::vector<double> &x, double &f, std::vector<double> lowerBounds,
            std::vector<double> upperBounds, std::vector<double> stepSize, double xTol, double fTol,
            int maxEval, GHSP_Data &_data);
void fullRun(const std::vector<double> x, GHSP_Data &data);
void outputVector(std::vector<unsigned int> const &set);
void outputVector(std::vector<float> const &set);

int main(int argc, char **argv) {
    // printf("Begin Index Construction and Search... \n");
    CLI::App app{"Index Construction"};
    std::string dataDirectory = "/users/cfoste18/scratch/sisap_data/";
    std::string dataset = "uniform";
    unsigned int dimension = 2;
    unsigned int datasetSize = 0;
    unsigned int testsetSize = 100;
    int numThreads = 1;
    bool verbose = true;
    app.add_option("-d,--dataset", dataset, "dataset");
    app.add_option("-D,--dimension", dimension, "Dimension");
    app.add_option("-N,--datasetSize", datasetSize, "Dataset Size");
    app.add_option("-T,--testsetSize", testsetSize, "Testset Size");
    app.add_option("-n,--numThreads", numThreads, "Number of Threads");
    app.add_option("-v,--verbose", verbose, "Show hRNG Statistics During Eval");

    // optimization parameters
    std::vector<double> radiusVector{};
    std::vector<double> stepSizeVector{};
    double xTol = 0.001;
    double fTol = 0.001;
    app.add_option("-r,--radiusVector", radiusVector, "starting radius vector")->required();
    app.add_option("-s,--stepSize", stepSizeVector, "starting radius vector")->required();
    app.add_option("-x,--xTol", xTol, "Stopping Condition: on radii values");
    app.add_option("-f,--fTol", fTol, "Stopping Condition: on function evaluation values");
    CLI11_PARSE(app, argc, argv);
    unsigned int n = radiusVector.size();

    // get dataset (with test set added)
    // dataset: 0...datsetSize-1, queryset: datasetSize-datasetSize+testsetSize
    // ALSO HAVE A VALIDATION SET
    unsigned int totalSize = datasetSize + testsetSize + testsetSize;
    float *dataPointer = NULL;
    if (dataset == "uniform") {  // validation set
        Datasets::uniformDataset(dataPointer, dimension, totalSize, 3);
    } else if (dataset == "LA") { 
        std::string data_path = "/users/cfoste18/scratch/LA/";
        Datasets::LA(data_path, dataPointer, dimension, totalSize);
        if (datasetSize == 0) {
            datasetSize = totalSize - 2*testsetSize;
        } else if (datasetSize + 2*testsetSize > totalSize) {
            datasetSize = totalSize - 2*testsetSize;
        }
    } else {
        printf("Unrecognized dataset: %s\n", dataset.c_str());
        return 0;
    }

    // Perform Optimization by COBYLA
    double fmin;
    std::vector<double> lowerBounds(n, 0);
    std::vector<double> upperBounds(n, (double) 100*radiusVector[0]);
    int maxEval = 200;
    GHSP_Data g_data(dimension, dataPointer, datasetSize, testsetSize);

    // Perform optimization
    COBLYA(n, radiusVector, fmin, lowerBounds, upperBounds, stepSizeVector, xTol, fTol, maxEval, g_data);

    // okay now do the real thing
    fullRun(radiusVector, g_data);

    printf("Done! Have a good day! \n");
    delete[] dataPointer;
    return 0;
}

/**
 * @brief
 *
 * @param n                 // num radii to optimize
 * @param x                 // starting values of the radii, returns optimal
 * @param f                 // returns optimal value
 * @param lowerBounds       // lower bound vector for each radii
 * @param upperBounds       // upper bound vector for each radii
 * @param stepSize          // starting step size vector for each radii
 * @param xTol              // stopping condition on radii vals
 * @param fTol              // stopping condition for eval
 * @param maxEval           // stoppint condition
 * @return true             // successful opt
 * @return false            // unsuccessful
 */
bool COBLYA(unsigned n, std::vector<double> &x, double &f, std::vector<double> lowerBounds,
            std::vector<double> upperBounds, std::vector<double> stepSize, double xTol, double fTol,
            int maxEval, GHSP_Data &_data) {
    printf("Using NLOPT Local Optimization Functions...\n");
    bool success = true;

    // convert hRNG_Data to void* for NLOPT
    void *ghsp_data = (void *)&_data;

    // create optimization object
    nlopt::opt opt = nlopt::opt(nlopt::LN_COBYLA, n);
    opt.set_lower_bounds(lowerBounds);  // sets the lower bounds
    opt.set_upper_bounds(upperBounds);  // sets the upper bounds
    opt.set_initial_step(stepSize);
    opt.set_min_objective(objectiveFunction, ghsp_data);
    opt.set_ftol_rel(fTol);
    opt.set_xtol_rel(xTol);
    opt.set_maxeval(maxEval);

    // perform optimization
    printf("Performing COBYLA Optimization...\n");
    printf("D, N, r, ER, |P|, Index Dist., Index Time (s), GHSP Distances, GHSP Time (ms)\n");
    std::chrono::high_resolution_clock::time_point tStart, tEnd;
    tStart = std::chrono::high_resolution_clock::now();
    int resultInt = opt.optimize(x, f);
    tEnd = std::chrono::high_resolution_clock::now();
    double optimizationTime = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();
    printf("End Optimization: Time=%.4f\n", optimizationTime);

    // print/save optimization termination reasoning
    std::string resultDescription = std::string("Terminated by Command ").append(std::to_string(resultInt));
    switch (resultInt) {
        case 1:
            // saver.printLine(resultDescription.append(": Success."));
            printf("%s: success!\n", resultDescription.c_str());
            break;
        case 2:
            // saver.printLine(resultDescription.append(": stopVal Reached."));
            printf("%s: stopVal Reached!\n", resultDescription.c_str());
            break;
        case 3:
            // saver.printLine(resultDescription.append(": fTol Reached."));
            printf("%s: fTol Reached!\n", resultDescription.c_str());
            break;
        case 4:
            // saver.printLine(resultDescription.append(": xTol Reached."));
            printf("%s: xTol Reached!\n", resultDescription.c_str());
            break;
        case 5:
            // saver.printLine(resultDescription.append(": maxEval Reached."));
            printf("%s: maxEval Reached!\n", resultDescription.c_str());
            break;
        case 6:
            // saver.printLine(resultDescription.append(": maxTime Reached."));
            printf("%s:  maxTime Reached!\n", resultDescription.c_str());
            break;
        case -1:
            // saver.printLine(resultDescription.append(": Generic Code Failure."));
            printf("%s: Generic Code Failure!\n", resultDescription.c_str());
            success = false;
            break;
        case -2:
            // saver.printLine(resultDescription.append(": Invalid Arguments."));
            printf("%s: Invalid Arguments!\n", resultDescription.c_str());
            success = false;
            break;
        case -3:
            // saver.printLine(resultDescription.append(": Ran out of Memory."));
            printf("%s: Ran out of Memory!\n", resultDescription.c_str());
            success = false;
            break;
        case -4:
            // saver.printLine(resultDescription.append(": Roundoff Errors."));
            printf("%s: Roundoff Errors!\n", resultDescription.c_str());
            success = false;
            break;
        case -5:
            // saver.printLine(resultDescription.append(": Forced Termination."));
            printf("%s: Forced Termination!\n", resultDescription.c_str());
            success = false;
            break;
        default:
            break;
    }

    return success;
}

/**
 * @brief Objective function to minimize for NLOPT. runs hRNG, looking to optimize search distance computations
 *
 * @param x
 * @param grad
 * @param hRNGData_Void
 * @return double
 */
double objectiveFunction(const std::vector<double> &x, std::vector<double> &grad, void *data_ghsp) {
    // convert void pointer back to hRNG_Data
    GHSP_Data *data = (GHSP_Data *)data_ghsp;
    unsigned int const dimension = data->_dimension;
    unsigned int const datasetSize = data->_datasetSize;
    unsigned int const testsetSize = data->_testsetSize;

    // get the radius vector
    std::vector<float> radiusVector(x.begin(), x.end());
    radiusVector.push_back(0.0f);
    int numLayers = radiusVector.size() + 1;

    // initialize sparse matrix for distance computations
    std::shared_ptr<SparseMatrix> sparseMatrix =
        std::make_shared<SparseMatrix>(data->_dataPointer, datasetSize + 2 * testsetSize, dimension);
    sparseMatrix->_datasetSize = datasetSize;

    //> Construct the pivot-index
    std::vector<PivotLayer> coverTree{};
    unsigned long long int dStart, dEnd;
    std::chrono::high_resolution_clock::time_point tStart, tEnd;
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    PivotIndex::CoverTree_Greedy(radiusVector, *sparseMatrix, coverTree);
    tEnd = std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    unsigned long long int distances_ps = (dEnd - dStart);
    double time_ps = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();

    //> Get Info on Hierarchy
    std::vector<unsigned int> pivotsPerLayer{};
    for (int i = 0; i < numLayers - 1; i++) {
        pivotsPerLayer.push_back((unsigned int)coverTree[i].pivotIndices->size());
    }
    pivotsPerLayer.push_back((unsigned int)datasetSize);

    //> Perform GHSP search on the pivot-index
    std::vector<std::vector<unsigned int>> neighbors{};
    neighbors.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::Hierarchical_HSP_Test(datasetSize + testsetSize + queryIndex, coverTree, *sparseMatrix, neighbors[queryIndex]);
    }
    tEnd = std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_hsp_pivot = (dEnd - dStart) / ((double)testsetSize);
    double time_hsp_pivot =
        std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double)testsetSize);

    printf("%u,%u,", dimension, datasetSize);
    outputVector(radiusVector);
    outputVector(pivotsPerLayer);
    printf("%llu,%.4f,", distances_ps, time_ps);
    printf("%.2f,%.4f\n", distances_hsp_pivot, time_hsp_pivot * 1000);

    return distances_hsp_pivot;
}

/**
 * @brief
 *
 * @param x             // optimal radii
 * @param data_ghsp     // data for everything
 */
void fullRun(const std::vector<double> x, GHSP_Data &data) {
    //> convert void pointer back to hRNG_Data
    unsigned int const dimension = data._dimension;
    unsigned int const datasetSize = data._datasetSize;
    unsigned int const testsetSize = data._testsetSize;

    //> get the radius vector
    std::vector<float> radiusVector(x.begin(), x.end());
    int numLayers = radiusVector.size() + 1;

    //> initialize sparse matrix for distance computations
    std::shared_ptr<SparseMatrix> sparseMatrix =
        std::make_shared<SparseMatrix>(data._dataPointer, datasetSize + 2 * testsetSize, dimension);
    sparseMatrix->_datasetSize = datasetSize;

    //> construct the pivot-index
    std::vector<PivotLayer> coverTree{};
    unsigned long long int dStart, dEnd;
    std::chrono::high_resolution_clock::time_point tStart, tEnd;
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    PivotIndex::CoverTree_Greedy(radiusVector, *sparseMatrix, coverTree);
    tEnd = std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    unsigned long long int distances_ps = (dEnd - dStart);
    double time_ps = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count();

    //> Get Info on Hierarchy
    std::vector<unsigned int> pivotsPerLayer{};
    for (int i = 0; i < numLayers - 1; i++) {
        pivotsPerLayer.push_back((unsigned int)coverTree[i].pivotIndices->size());
    }
    pivotsPerLayer.push_back((unsigned int)datasetSize);

    // perform GHSP search on the pivot-index
    std::vector<std::vector<unsigned int>> results_hsp_pivot{};
    results_hsp_pivot.resize(testsetSize);
    dStart = sparseMatrix->_distanceComputationCount;
    tStart = std::chrono::high_resolution_clock::now();
    for (unsigned int queryIndex = 0; queryIndex < testsetSize; queryIndex++) {
        GHSP::Hierarchical_HSP_Test(datasetSize + queryIndex, coverTree, *sparseMatrix, results_hsp_pivot[queryIndex]);
    }
    tEnd = std::chrono::high_resolution_clock::now();
    dEnd = sparseMatrix->_distanceComputationCount;
    double distances_hsp_pivot = (dEnd - dStart) / ((double)testsetSize);
    double time_hsp_pivot =
        std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart).count() / ((double)testsetSize);


    printf("----------------------------------------------\n");
    printf("N,R,|P|,Index Distances,Index Time (s), Pivot HSP Distances, Pivot HSP Time (ms) \n");
    printf("%u,", datasetSize);
    outputVector(radiusVector);
    outputVector(pivotsPerLayer);
    printf("%llu,%.4f,", distances_ps, time_ps);
    printf("%.2f,%.4f,", distances_hsp_pivot, time_hsp_pivot * 1000);
    printf("\n");

    return;
}

void outputVector(std::vector<unsigned int> const &set) {
    printf("{;");
    std::vector<unsigned int>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%u;", (*it1));
    }
    printf("},");
}

void outputVector(std::vector<float> const &set) {
    printf("{;");
    std::vector<float>::const_iterator it1;
    for (it1 = set.begin(); it1 != set.end(); it1++) {
        printf("%.6f;", (*it1));
    }
    printf("},");
}
