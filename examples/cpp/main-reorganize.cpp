/*
    Cole Foster
    2023-07-21
    Hierarchical HSP
*/
#include <chrono>
#include <cstdio>
#include <vector>

#include "CLI11.hpp"
#include "datasets.hpp"

#include "HierarchicalHSP.hpp"

void print_vector(std::vector<elementID> const& vec);

int main(int argc, char **argv) {
    CLI::App app{""};
    size_t dimension = 2;
    size_t dataset_size = 1000;
    app.add_option("-D,--dimension", dimension, "dimension");
    app.add_option("-N,--dataset_size", dataset_size, "dataset size");
    CLI11_PARSE(app, argc, argv);

    // retrieve the dataset
    float* data_pointer;
    Datasets::uniformDataset(data_pointer, dimension, dataset_size, 3);
    printf("d=uniform,N=%u,D=%u\n", dataset_size, dimension);

    // initialize the hhsp class
    HierarchicalHSP HHSP(dimension, dataset_size, data_pointer);

    //> Find Some HSP Neighbors
    elementID index1 = 10;
    std::vector<elementID> neighbors_index1{};
    HHSP.HSP_Test(index1, neighbors_index1);
    printf("HSP(%u) = ",index1);
    print_vector(neighbors_index1);

    elementID index2 = 14;
    std::vector<elementID> neighbors_index2{};
    HHSP.HSP_Test(index2, neighbors_index2);
    printf("HSP(%u) = ",index2);
    print_vector(neighbors_index2);

    // // distance comp
    // elementID index2 = 14;
    // HHSP.print_representation(index1);
    // HHSP.print_representation(index2);
    // float const distance = HHSP.computeDistance(index1,index2);
    // printf("Euclidean Distance: %.6f\n",distance);


    printf("Done! Have a good day! \n");
    delete[] data_pointer;
    return 0;
}

void print_vector(std::vector<elementID> const& vec) {
    printf("{");
    for (int i = 0; i < (int) vec.size(); i++) {
        printf("%u,",(unsigned int) vec[i]);
    }
    printf("}\n");
    return;
}