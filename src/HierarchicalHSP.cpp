#include "HierarchicalHSP.hpp"
#include <queue>

void HierarchicalHSP::HSP_Test(float* query_pointer, std::vector<elementID>& neighbors) const {
    neighbors.clear();

    // only perform on k closest elements
    std::vector<std::pair<float, elementID>> active_list{};
    active_list.reserve(dataset_size_);

    // find next nearest neighbor and create list of distances
    elementID index1;
    float distance_Q1 = HUGE_VAL;
    for (elementID index = 0; index < (elementID)dataset_size_; index++) {
        float const distance = compute_distance(query_pointer, internalDataPointer_(index));
        if (distance <= 0) {
            neighbors.push_back(index);
            continue;
        }
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
        for (int it2 = 0; it2 < (int)active_list_copy.size(); it2++) {
            elementID const index2 = active_list_copy[it2].second;
            float const distance_Q2 = active_list_copy[it2].first;
            if (index2 == index1) continue;
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
//                     fstdistfunc(dataPointer + queryIndex * dimension, dataPointer + index2 * dimension,
//                     void_dimension);
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

void HierarchicalHSP::create_index(std::vector<float> radius_vector) {
    radius_vector.push_back(0);  // bottom layer has zero-radius
    int const num_layers = (int)radius_vector.size();

    // //> Change to Effective Radius Vector -> Makes them Covering Radii for entire pivot domains
    //          ** this should be done ahead of time, not here
    // std::vector<float> temp_vec = radiusVector;
    // for (int i = 0; i < (int)radiusVector.size(); i++) {
    //     for (int j = i + 1; j < (int)temp_vec.size(); j++) {
    //         radiusVector[i] += (float)temp_vec[j];
    //     }
    // }

    //> Initialize Each Layer of the Cover Tree
    (*index_).clear();
    (*index_).resize(num_layers - 1);
    for (int layer_index = 0; layer_index < num_layers - 1; layer_index++) {
        (*index_)[layer_index] = PivotLayer(layer_index, num_layers, radius_vector[layer_index]);
    }

    //> Incremental construction in a top-down fashion
    for (elementID query_index = 0; query_index < (elementID)dataset_size_; query_index++) {
        int orphan_layer = 0;        // 1 Below Lowest Layer With a Parent
        unsigned int lowest_parent;  // Closest Parent in the Layer Above Orphanage

        //> Top-Down approach, find the lowest layer with a parent
        std::vector<unsigned int> grandparent_list{};
        for (int layer_index = 0; layer_index <= num_layers - 1; layer_index++) {
            //> Collect the Spotlight |--> Pivots in layer local to Q (may be a parent)
            std::vector<unsigned int> spotlight{};
            if (layer_index == 0) {
                // Collect All Top Layer Pivots
                tsl::sparse_set<unsigned int> const& top_layer_pivots = *((*index_)[0].get_pivotIndices_ptr());
                spotlight.insert(spotlight.end(), top_layer_pivots.begin(), top_layer_pivots.end());

            } else {
                // Collect the Children of "Grandparents", may be parents of Q
                for (int it1 = 0; it1 < (int)grandparent_list.size(); it1++) {
                    unsigned int const grandparent_pivot = grandparent_list[it1];
                    std::vector<unsigned int> const& pivot_domain =
                        (*index_)[layer_index - 1].get_pivotChildren(grandparent_pivot);
                    spotlight.insert(spotlight.end(), pivot_domain.begin(), pivot_domain.end());
                }
            }

            //> Find the closest parent of Q from the Spotlight
            float closest_parent_distance = HUGE_VAL;
            grandparent_list.clear();  // grandparents: any pivot within the covering radius r of Q
            for (int it2 = 0; it2 < (int)spotlight.size(); it2++) {
                unsigned int const candidate_parent_pivot = spotlight[it2];
                float const distance = compute_distance(query_index, candidate_parent_pivot);

                // Test if Q is within the Covering Radius
                if (distance <= radius_vector[layer_index]) {
                    grandparent_list.push_back(candidate_parent_pivot);  // could have children as lower parent

                    // Test if Q is an actual parent of Q
                    if (distance <= radius_vector[layer_index] - radius_vector[layer_index + 1]) {
                        // is this new lowest parent layer?
                        if (orphan_layer < layer_index + 1) {
                            orphan_layer = layer_index + 1;  // no parent yet in layer below
                        }

                        // is this the closest parent?
                        if (distance < closest_parent_distance) {
                            closest_parent_distance = distance;
                            lowest_parent = candidate_parent_pivot;
                        }
                    }
                }
            }

            //> Break if no potential future parents
            if (grandparent_list.empty()) break;
        }

        //> Nesting: Add Query To Orphan Layer and All Layers Below
        for (int layer_index = orphan_layer; layer_index < num_layers; layer_index++) {
            //
            //> Becomes a Pivot on Every Lower Layer Except Bottom: Not Explicitly Represented
            if (layer_index < num_layers - 1) {
                (*index_)[layer_index].addPivot(query_index);
            }

            //> Becomes a Child Of Itself
            if (layer_index > 0) {
                (*index_)[layer_index - 1].addChild(query_index, query_index);
            }

            //> Becomes a Parent Of Itself
            if (layer_index < num_layers - 1 && layer_index > orphan_layer) {
                (*index_)[layer_index].addParent(query_index, query_index);
            }
        }

        // Add The Parent of the Query, Recursively Update Parents MaxChildDistance
        if (orphan_layer > 0) {
            if (orphan_layer < num_layers - 1) {
                (*index_)[orphan_layer].addParent(query_index, lowest_parent);
            }
            (*index_)[orphan_layer - 1].addChild(lowest_parent, query_index);
            (*index_)[orphan_layer - 1].updateMaxChildDistance(lowest_parent,
                                                               compute_distance(lowest_parent, query_index));

            // Update Max Child Distance for All Parents of Parents
            unsigned int parent = lowest_parent;
            for (int layer_index = orphan_layer - 1; layer_index > 0; layer_index--) {
                unsigned int grandParent = (*index_)[layer_index].get_parentIndex(parent);
                (*index_)[layer_index - 1].updateMaxChildDistance(grandParent,
                                                                  compute_distance(grandParent, query_index));
                parent = grandParent;
            }
        }
    }

    return;
};

/**
 * @brief Validate the Coverage and Nesting Properties of the Cover Tree, but not Separation.
 *
 * @return true
 * @return false
 */
bool HierarchicalHSP::validate_index() {
    int const num_layers = (*index_).size() + 1;

    // print stats:
    printf("LayerID, Radius, NumPivots\n");
    for (int layer_index = 0; layer_index < num_layers - 1; layer_index++) {
        printf("  %d: %.4f, %u\n", layer_index, (*index_)[layer_index].radius, (*index_)[layer_index].get_num_pivots());
    }
    printf("  %d: %.4f, %u\n", num_layers - 1, 0.0000, dataset_size_);

    // each point should be represented
    std::vector<elementID> points_list{};
    points_list.resize(dataset_size_);
    for (elementID query_index = 0; query_index < (elementID)dataset_size_; query_index++) {
        points_list[query_index] = query_index;
    }

    // Top-Layer Iteration
    tsl::sparse_set<unsigned int> const& pivotIndices = *((*index_)[0].get_pivotIndices_ptr());
    tsl::sparse_set<unsigned int>::const_iterator it1;
    for (it1 = pivotIndices.begin(); it1 != pivotIndices.end(); it1++) {
        unsigned int const pivot_index = (*it1);
        std::vector<unsigned int> parents_list{};
        parents_list.resize(num_layers);
        parents_list[0] = pivot_index;
        bool flag_success = recursive_depth_first_check(parents_list, 0, points_list);
        if (!flag_success) return false;
    }

    // should be no points left
    if (points_list.size() > 0) {
        printf("There are points not covered by the index! Coverage violated!\n");
        printf("  -> %u: ", points_list.size());
        print_vector(points_list);
        return false;
    }

    printf("Index Ensures Minimal Coverage of the Dataset!\n");
    return true;
}

/**
 * @brief Recursive Function for Performing Depth-First Cover Tree Validation
 *
 * @param parents
 * @param layer_index
 * @param points_list
 * @return true
 * @return false
 */
bool HierarchicalHSP::recursive_depth_first_check(std::vector<unsigned int> parents, int const layer_index,
                                                  std::vector<unsigned int>& points_list) {
    int num_layers = (*index_).size() + 1;

    // Iterate through the Domain of Lowest Parent
    std::vector<unsigned int> const& pivot_domain = (*index_)[layer_index].get_pivotChildren(parents[layer_index]);
    for (int it1 = 0; it1 < (int)pivot_domain.size(); it1++) {
        unsigned int const child_index = pivot_domain[it1];
        float child_radius = 0;        // radius of the child pivot
        float child_max_distance = 0;  // max distance to child of child
        if (layer_index + 1 < num_layers - 1) {
            child_radius = (*index_)[layer_index + 1].radius;
            // child_max_distance = (*index_)[layer_index + 1].get_maxChildDistance(child_index);
        }

        // check child_index within radius/maxDistance of every parent index
        for (int coarser_layer_index = 0; coarser_layer_index <= layer_index; coarser_layer_index++) {
            unsigned int const coarser_parent = parents[coarser_layer_index];
            float const pivot_radius = (*index_)[coarser_layer_index].radius;
            float const pivot_max_child_distance = (*index_)[coarser_layer_index].get_maxChildDistance(coarser_parent);

            float const distance = compute_distance(coarser_parent, child_index);
            if (distance > pivot_radius - child_radius) {
                printf("child_index %u on L-%d not within radius of ParentIndex %u on L-%d!\n", child_index,
                       layer_index + 1, coarser_parent, coarser_layer_index);
                return false;
            }

            if (distance > pivot_max_child_distance) {
                printf(
                    "child_index %u on L-%d -max child distance- not within max-child-distance of ParentIndex %u on "
                    "L-%d!\n",
                    child_index, layer_index + 1, coarser_parent, coarser_layer_index);
                return false;
            }
        }

        if (layer_index + 1 < num_layers - 1) {
            parents[layer_index + 1] = child_index;
            recursive_depth_first_check(parents, layer_index + 1, points_list);
        } else {
            // if bottom layer, remove points from list
            std::vector<unsigned int>::iterator it3a = std::find(points_list.begin(), points_list.end(), child_index);
            if (it3a == points_list.end()) {
                printf("Child has already been represented! Coverage not minimal!\n");
                printf("c=%u\n", child_index);
                return false;
            } else {
                points_list.erase(it3a);
            }
        }
    }

    return true;
}

// for ranked lists
typedef std::priority_queue<std::pair<float, elementID>, std::vector<std::pair<float, elementID>>,
                            std::greater<std::pair<float, elementID>>>
    PriorityQueue;

/**
 * @brief Perform the HSP Test Hierarchically Given a Multi-Layer Cover Tree
 *
 * @param query_pointer
 * @param neighbors
 */
void HierarchicalHSP::Hierarchical_HSP_Test(float* query_pointer, std::vector<elementID>& neighbors) const {
    neighbors.clear();
    int const num_layers = (*index_).size() + 1;

    //  - initialize query distance storage
    elementID query_index = dataset_size_ + 1;  // FAKE INDEX: for distance storages
    tsl::sparse_map<int, std::vector<float>> distance_store_array{};
    distance_store_array[query_index] = std::vector<float>(dataset_size_, -1);

    //  - initialize Lists for Multi-Layer HSP Search
    std::vector<PriorityQueue> active_pivots{};  // pivots safe from invalidation by any existing neighbors
    active_pivots.resize(num_layers);
    std::vector<PriorityQueue> intermediate_pivots{};  // pivots questionable about safety
    intermediate_pivots.resize(num_layers);

    //  - compute distance to all top layer pivots to initialize the queue
    {
        tsl::sparse_set<unsigned int> const& pivot_indices = *((*index_)[0].get_pivotIndices_ptr());
        tsl::sparse_set<unsigned int>::const_iterator it1;
        for (it1 = pivot_indices.begin(); it1 != pivot_indices.end(); it1++) {
            elementID const pivot_index = (elementID)(*it1);
            float const distance = get_distance(query_index, query_pointer, pivot_index, distance_store_array);
            active_pivots[0].push(std::make_pair(distance, pivot_index));
        }
    }

    //  - count the total number of points in lists
    unsigned int num_remaining_points = 0;
    for (int layer_index = 0; layer_index < num_layers; layer_index++) {
        num_remaining_points += (active_pivots[layer_index].size() + intermediate_pivots[layer_index].size());
    }

    // Begin the Cover Tree- HSP Algorithm Loop
    while (num_remaining_points >= 0) {
        /**
         * ======================================================
         *   STEP 1: Active Nearest Neighbor Search
         *     - find index1 as the next hsp neighbor of Q
         *
         * ======================================================
         */
        elementID index1;
        float dmin = HUGE_VAL;

        //> Check the active pivots bottom-up
        for (int layer_index = num_layers - 1; layer_index >= 0; layer_index--) {
            PriorityQueue list_copy = active_pivots[layer_index];

            // iterate through this copy of the queue
            std::pair<float, elementID> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                elementID const pivot_index = topOfQueue.second;
                float const distance = topOfQueue.first;
                list_copy.pop();

                if (distance < dmin) {
                    if (std::find(neighbors.begin(), neighbors.end(), pivot_index) != neighbors.end()) continue;
                    dmin = distance;
                    index1 = pivot_index;
                } else {
                    break;
                }
            }
        }

        //> Check the intermediate pivots bottom-up
        for (int layer_index = num_layers - 1; layer_index >= 0; layer_index--) {
            PriorityQueue list_copy = intermediate_pivots[layer_index];

            // iterate through this copy of the queue
            std::pair<float, elementID> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                elementID const pivot_index = topOfQueue.second;
                float const distance = topOfQueue.first;
                list_copy.pop();

                if (distance < dmin) {
                    if (std::find(neighbors.begin(), neighbors.end(), pivot_index) != neighbors.end()) continue;

                    // update dmin if p2 is valid
                    bool flag_validated =
                        validate_point(pivot_index, query_index, query_pointer, neighbors, distance_store_array);
                    if (flag_validated) {
                        dmin = distance;
                        index1 = pivot_index;
                    }
                } else {
                    break;
                }
            }
        }

        //> Check the Domains of Active Pivots Top-Down
        std::vector<PriorityQueue> copy_active_pivots = active_pivots;
        for (int layer_index = 0; layer_index < num_layers - 1; layer_index++) {
            //
            // iterate through the list of pivots on this layer
            std::pair<float, elementID> topOfQueue;
            while (copy_active_pivots[layer_index].size() > 0) {
                topOfQueue = copy_active_pivots[layer_index].top();
                copy_active_pivots[layer_index].pop();

                // get pivot information
                elementID const pivot_index = topOfQueue.second;
                float const distance_Qp = topOfQueue.first;
                float const radius = (*index_)[layer_index].get_maxChildDistance(pivot_index);

                // check if domain may contain a closer point
                if (distance_Qp <= dmin + radius) {
                    //
                    // iterate through the pivot domain
                    std::vector<unsigned int> const& pivot_domain =
                        (*index_)[layer_index].get_pivotChildren(pivot_index);
                    for (int it1 = 0; it1 < pivot_domain.size(); it1++) {
                        elementID const child_index = (elementID)pivot_domain[it1];
                        float const distance_Qc =
                            get_distance(query_index, query_pointer, child_index, distance_store_array);

                        // check if this child is the new closest point
                        if (distance_Qc < dmin) {
                            if (std::find(neighbors.begin(), neighbors.end(), child_index) == neighbors.end()) {
                                dmin = distance_Qc;
                                index1 = child_index;
                            }
                        }

                        // check if domain of child may contain a closer point
                        if (layer_index < num_layers - 2) {
                            float const childRadius = (*index_)[layer_index + 1].get_maxChildDistance(child_index);
                            if (distance_Qc <= dmin + childRadius) {
                                copy_active_pivots[layer_index + 1].push(std::make_pair(distance_Qc, child_index));
                            }
                        }
                    }
                }
            }
        }

        //> Check the Domains of Intermediate Pivots Top-Down
        //  - improvement to make: if we have to examine the entire domain, then remove parent pivot from list
        //  - improvement to make: we can check if the entire domain is safe/invalid, then move to active list
        for (int layer_index = 0; layer_index < num_layers - 1; layer_index++) {
            PriorityQueue list_copy = intermediate_pivots[layer_index];
            intermediate_pivots[layer_index] = PriorityQueue();

            // iterate through the list of pivots on this layer
            std::pair<float, elementID> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                list_copy.pop();

                // get pivot information
                elementID const pivot_index = topOfQueue.second;
                float const distance_Qp = topOfQueue.first;
                float const radius = (*index_)[layer_index].get_maxChildDistance(pivot_index);

                // check if domain may contain a closer point
                if (distance_Qp <= dmin + radius) {
                    //
                    // iterate through the pivot domain
                    std::vector<unsigned int> const& pivot_domain =
                        (*index_)[layer_index].get_pivotChildren(pivot_index);
                    for (int it1 = 0; it1 < pivot_domain.size(); it1++) {
                        elementID const child_index = (elementID)pivot_domain[it1];
                        float const distance_Qc =
                            get_distance(query_index, query_pointer, child_index, distance_store_array);

                        // if bottom layer, no domain
                        if (layer_index + 1 == num_layers - 1) {
                            if (std::find(neighbors.begin(), neighbors.end(), child_index) != neighbors.end()) continue;

                            // if point is valid, it is active (not intermediate) and can update dmin
                            bool flag_validated = validate_point(child_index, query_index, query_pointer, neighbors,
                                                                 distance_store_array);
                            if (flag_validated) {
                                if (distance_Qc < dmin) {
                                    dmin = distance_Qc;
                                    index1 = child_index;
                                }
                                active_pivots[layer_index + 1].push(std::make_pair(distance_Qc, child_index));
                            }
                        } else {
                            float const child_radius = (*index_)[layer_index + 1].get_maxChildDistance(child_index);

                            // validate entire pivot domain
                            //
                            //  -----------------
                            //
                            int pivot_outcome = 2;
                            bool flag_validated = validate_pivot(child_index, child_radius, query_index, query_pointer,
                                                                 neighbors, distance_store_array, pivot_outcome);
                            if (flag_validated) {
                                if (distance_Qc < dmin) {
                                    if (std::find(neighbors.begin(), neighbors.end(), child_index) == neighbors.end()) {
                                        dmin = distance_Qc;
                                        index1 = child_index;
                                    }
                                }
                            }

                            // based on entire pivot domain validation
                            if (pivot_outcome == 0) {                      // active, not intermediate pivot
                                if (distance_Qc <= dmin + child_radius) {  // must investigate lower domain
                                    intermediate_pivots[layer_index + 1].push(std::make_pair(distance_Qc, child_index));
                                } else {
                                    active_pivots[layer_index + 1].push(std::make_pair(distance_Qc, child_index));
                                }
                            } else if (pivot_outcome == 1) {  // pivot and domain eliminated
                                continue;

                            } else if (pivot_outcome == 2) {  // intermediate pivot and domain
                                intermediate_pivots[layer_index + 1].push(std::make_pair(distance_Qc, child_index));
                            }
                        }
                    }
                } else {
                    intermediate_pivots[layer_index].push(std::make_pair(distance_Qp, pivot_index));
                }
            }
        }

        // assign the new hsp neighbor
        if (dmin > 10000) break;  // large, magic number signaling no new neighbor possible
        elementID const x1 = index1;
        float const distance_Q1 = get_distance(query_index, query_pointer, x1, distance_store_array);
        float* x1_pointer = internalDataPointer_(x1);
        neighbors.push_back(x1);
        distance_store_array[x1] = std::vector<float>(dataset_size_, -1);  // store dist. for neighbors

        /**
         * =================================================================
         *   STEP 2: Validation of Active Points
         *     - remove points that are invalidated by the new HSP neighbor
         *
         * =================================================================
         */

        //> Iterate through the intermediate list in each layer, top-down.
        for (int layer_index = 0; layer_index < num_layers; layer_index++) {
            // start by iterating through the intermediate pivots
            PriorityQueue list_copy = intermediate_pivots[layer_index];
            intermediate_pivots[layer_index] = PriorityQueue();
            std::pair<float, elementID> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                list_copy.pop();

                // get pivot information
                elementID const pivot_index = topOfQueue.second;
                float const distance_Q2 = topOfQueue.first;
                float radius = 0.0f;
                if (layer_index < num_layers - 1) {
                    radius = (*index_)[layer_index].get_maxChildDistance(pivot_index);
                }

                //> Pivot sufficiently close to Q. Add domain to lower layer.
                if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                    //  - add the entire domain to the intermeidate list of the lower layer
                    if (layer_index < num_layers - 1) {
                        std::vector<unsigned int> const& pivot_domain =
                            (*index_)[layer_index].get_pivotChildren(pivot_index);
                        for (int it1 = 0; it1 < pivot_domain.size(); it1++) {
                            elementID const child_index = (elementID)pivot_domain[it1];
                            float const distance_Qc =
                                get_distance(query_index, query_pointer, child_index, distance_store_array);
                            intermediate_pivots[layer_index + 1].push(std::make_pair(distance_Qc, child_index));
                        }
                    }

                    //  - add the pivot back to the list
                    else {
                        intermediate_pivots[layer_index].push(std::make_pair(distance_Q2, pivot_index));
                    }
                }

                //> Check pivot for invalidation
                else {
                    float const distance_12 = get_distance(x1, x1_pointer, pivot_index, distance_store_array);

                    //  - proposition: entire domain invalidated
                    if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                        continue;
                    }
                    //  - add the pivot back to the list
                    else {
                        intermediate_pivots[layer_index].push(std::make_pair(distance_Q2, pivot_index));
                    }
                }
            }
        }

        //> Iterate through the intermediate list in each layer, top-down.
        for (int layer_index = 0; layer_index < num_layers; layer_index++) {
            // start by iterating through the intermediate pivots
            PriorityQueue list_copy = active_pivots[layer_index];
            active_pivots[layer_index] = PriorityQueue();
            std::pair<float, elementID> topOfQueue;
            while (list_copy.size() > 0) {
                topOfQueue = list_copy.top();
                list_copy.pop();

                // get pivot information
                elementID const pivot_index = topOfQueue.second;
                float const distance_Q2 = topOfQueue.first;
                float radius = 0;
                if (layer_index < num_layers - 1) {
                    radius = (*index_)[layer_index].get_maxChildDistance(pivot_index);
                }

                //> Pivot sufficiently close to Q. Add domain to lower layer list.
                if (distance_Q1 >= distance_Q2 - radius || distance_Q2 < radius) {
                    //  - add the entire domain to the intermeidate list of the lower layer
                    if (layer_index < num_layers - 1) {
                        std::vector<unsigned int> const& pivot_domain =
                            (*index_)[layer_index].get_pivotChildren(pivot_index);
                        for (int it1 = 0; it1 < pivot_domain.size(); it1++) {
                            elementID const child_index = (elementID)pivot_domain[it1];
                            float const distance_Qc =
                                get_distance(query_index, query_pointer, child_index, distance_store_array);

                            intermediate_pivots[layer_index + 1].push(std::make_pair(distance_Qc, child_index));
                        }
                    }

                    //  - add the pivot back to the list
                    else {
                        intermediate_pivots[layer_index].push(std::make_pair(distance_Q2, pivot_index));
                    }
                }

                //> Check pivot for invalidation
                else {
                    float const distance_12 = get_distance(x1, x1_pointer, pivot_index, distance_store_array);

                    //  - proposition: entire domain invalidated
                    if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                        continue;
                    }
                    //  - proposition: Entire Domain Safe
                    else if (distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                        active_pivots[layer_index].push(std::make_pair(distance_Q2, pivot_index));
                    }
                    //  - add the pivot back to the list
                    else {
                        intermediate_pivots[layer_index].push(std::make_pair(distance_Q2, pivot_index));
                    }
                }
            }
        }

        // update the number of remaining points
        for (int layer_index = 0; layer_index < num_layers; layer_index++) {
            num_remaining_points += (active_pivots[layer_index].size() + intermediate_pivots[layer_index].size());
        }
    }

    return;
}

// index1 should be the query or a neighbor. index2 always member of dataset
float HierarchicalHSP::get_distance(elementID const index1, float* index1_ptr, elementID const index2,
                                    tsl::sparse_map<int, std::vector<float>>& distance_store_array) const {
    float distance = distance_store_array.at(index1)[index2];
    if (distance < 0) {
        distance = compute_distance(index1_ptr, internalDataPointer_(index2));
        distance_store_array.at(index1)[index2] = distance;
    }
    return distance;
}

/**
 * @brief Check if the element point_index (in bottom layer) is invalidated by any existing neighbor
 *
 * @param point_index
 * @param query_index
 * @param query_ptr
 * @param neighbors
 * @param distance_store_array
 * @return true
 * @return false
 */
bool HierarchicalHSP::validate_point(elementID const point_index, elementID const query_index, float* query_ptr,
                                     std::vector<elementID> const& neighbors,
                                     tsl::sparse_map<int, std::vector<float>>& distance_store_array) const {
    // point_valid: flag
    //  - true: index of pivot is not invalidated by any existing HSP neighbor
    //  - false: it is invalidated
    bool point_valid = true;

    // HSP Test Against HSP Neighbors
    float const distance_Q2 = get_distance(query_index, query_ptr, point_index, distance_store_array);
    for (int it1 = 0; it1 < (int)neighbors.size(); it1++) {
        elementID const x1_bar = neighbors[it1];
        float const distance_Q1 = get_distance(query_index, query_ptr, x1_bar, distance_store_array);
        float const distance_12 = get_distance(x1_bar, internalDataPointer_(x1_bar), point_index, distance_store_array);

        // check point against normal HSP inequalities
        if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
            point_valid = false;
            break;
        }
    }

    return point_valid;
}

bool HierarchicalHSP::validate_pivot(elementID const pivot_index, float const radius, elementID const query_index,
                                     float* query_ptr, std::vector<elementID> const& neighbors,
                                     tsl::sparse_map<int, std::vector<float>>& distance_store_array,
                                     int& pivot_outcome) const {
    // point_valid: flag
    //  - true: index of pivot (on bottom layer) is not invalidated by any existing HSP neighbor
    //  - false: it is invalidated
    bool point_valid = true;

    // pivot_outcome:
    //      0: pivot domain fully safe by ALL neighbors, not intermediate
    //      1: pivot domain fully invalidated by a neighbor, remove from consideration
    //      2: neither, remains intermediate
    pivot_outcome = 0;  // default safe, interference/intermediate is found

    // iterate through each HSP neighbor, check against HHSP inequalities
    float const distance_Q2 = get_distance(query_index, query_ptr, pivot_index, distance_store_array);
    for (int it1 = 0; it1 < (int)neighbors.size(); it1++) {
        elementID const x1_bar = neighbors[it1];
        float const distance_Q1 = get_distance(query_index, query_ptr, x1_bar, distance_store_array);
        float const distance_12 = get_distance(x1_bar, internalDataPointer_(x1_bar), pivot_index, distance_store_array);

        // check point against normal HSP inequalities
        if (point_valid) {  // if not yet invalidated
            if (distance_Q1 < distance_Q2 && distance_12 < distance_Q2) {
                point_valid = false;
            }
        }

        // check if the HSP neighbor x1 fully invalidates the domain
        if (distance_Q1 < distance_Q2 - radius) {
            if (distance_Q2 * distance_Q2 - distance_12 * distance_12 > 2 * radius * distance_Q1) {
                point_valid = false;
                pivot_outcome = 1;
                break;
            }
        }

        // check if the domain is safe
        if (pivot_outcome == 0) {  // only check if if not yet intermediate
            if (distance_Q1 > distance_Q2 + radius ||
                distance_12 * distance_12 - distance_Q2 * distance_Q2 > 2 * radius * distance_Q1) {
                continue;  // safe, remains 0
            } else {
                pivot_outcome = 2;  // intermediate
            }
        }
    }

    return point_valid;
}