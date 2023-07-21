#include "HierarchicalHSP.hpp"


/**
 * @brief Given a neighborhood list L, find the HSP neighbors of Q
 *
 * 
 * improvements:
 *  - multithreading
 *  - when query is member of dataset, don't check for interference.
 *
 */
void HierarchicalHSP::HSP_Test(float* query_pointer, std::vector<elementID>& neighbors) const {
    neighbors.clear();

    // only perform on k closest elements
    std::vector<std::pair<float, elementID>> active_list{};
    active_list.resize(dataset_size_);

    // find next nearest neighbor and create list of distances
    elementID index1;
    float distance_Q1 = HUGE_VAL;
    for (elementID index = 0; index < (elementID) dataset_size_; index++) {
        float const distance = compute_distance(query_pointer, internalDataPointer_(index));
        if (distance < distance_Q1) {
            distance_Q1 = distance;
            index1 = index;
        }
        active_list.emplace_back(distance, index);
    }

    // perform the hsp loop
    while (active_list.size() > 0) {
        neighbors.push_back(index1);
        float* index1_ptr = internalDataPointer_(index1);

        //  - set up for the next hsp neighbor
        elementID index1_next;
        float distance_Q1_next = HUGE_VAL;

        //  - initialize the active_list for next iteration
        std::vector<std::pair<float, elementID>> active_list_copy = active_list;
        active_list.clear();

        //  - check each point for elimination
        for (int it1 = 0; it1 < (int)active_list_copy.size(); it1++) {
            elementID const index2 = active_list_copy[it1].second;
            if (index2 == index1) continue;
            float const distance_Q2 = active_list_copy[it1].first;
            float const distance_12 = compute_distance(index1_ptr, internalDataPointer_(index2));

            // check the hsp inequalities: add if not satisfied
            if (distance_Q1 >= distance_Q2 || distance_12 >= distance_Q2) {
                active_list.emplace_back(distance_Q2, index2);
                if (distance_Q2 < distance_Q1_next) {
                    distance_Q1_next = distance_Q2;
                    index1_next = index2;
                }
            }
        }

        // setup the next hsp neighbor
        index1 = index1_next;
        distance_Q1 = distance_Q1_next;
    }



    return;
}


// /**
//  * @brief Given a neighborhood list L, find the HSP neighbors of Q
//  *

//  */
// void HierarchicalHSP::HSP_Test(float* query_pointer, std::vector<elementID>& neighbors, bool multithreaded);
//     int const num_threads_to_use = 1;
//     if (multithreaded) num_threads_to_use = num_threads_;
//     neighbors.clear();

//     // find next nearest neighbor
//     unsigned int index1;
//     float distance_Q1 = HUGE_VAL;

// #pragma omp parallel num_threads(numThreads)
//     {
//         unsigned int index1a;
//         float distance_Q1a = HUGE_VAL;

// // compute distances in parallel
// #pragma omp for schedule(static)
//         for (unsigned int i1 = 0; i1 < L.size(); i1++) {
//             unsigned int index = L[i1];
//             if (index == queryIndex) continue;
//             float const distance =
//                 fstdistfunc(dataPointer + queryIndex * dimension, dataPointer + index * dimension, void_dimension);
//             if (distance < distance_Q1a) {
//                 distance_Q1a = distance;
//                 index1a = index;
//             }
//         }

// // update for nearest neighbor
// #pragma omp critical(updateNN)
//         {
//             if (distance_Q1a < distance_Q1) {
//                 distance_Q1 = distance_Q1a;
//                 index1 = index1a;
//             }
//         }
//     }

//     // now, eliminate points and find next hsp neighbors
//     while (L.size() > 0) {
//         neighbors.push_back(index1);
//         std::vector<unsigned int> L_copy = L;
//         L.clear();
//         float distance_Q1_temp = HUGE_VAL;

// // in parallel, eliminate points and find next neighbor
// #pragma omp parallel num_threads(numThreads)
//         {
//             unsigned int index1a;
//             float distance_Q1a = HUGE_VAL;
//             std::vector<unsigned int> La{};

// // compute distances in parallel
// #pragma omp for schedule(static)
//             for (unsigned int i1 = 0; i1 < L_copy.size(); i1++) {
//                 unsigned int const index2 = L_copy[i1];
//                 if (index2 == index1 || index2 == queryIndex) continue;
//                 float const distance_Q2 =
//                     fstdistfunc(dataPointer + queryIndex * dimension, dataPointer + index2 * dimension, void_dimension);
//                 float const distance_12 =
//                     fstdistfunc(dataPointer + index1 * dimension, dataPointer + index2 * dimension, void_dimension);

//                 // check inequalities
//                 if (distance_Q1 >= distance_Q2 || distance_12 >= distance_Q2) {
//                     La.push_back(index2);
//                     if (distance_Q2 < distance_Q1a) {
//                         distance_Q1a = distance_Q2;
//                         index1a = index2;
//                     }
//                 }
//             }

// // update for nearest neighbor
// #pragma omp critical(updateNN)
//             {
//                 L.insert(L.end(), La.begin(), La.end());
//                 if (distance_Q1a < distance_Q1_temp) {
//                     distance_Q1_temp = distance_Q1a;
//                     index1 = index1a;
//                 }
//             }
//         }
//         distance_Q1 = distance_Q1_temp;
//     }

//     return;
// }