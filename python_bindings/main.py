import sys
import argparse
import numpy as np
import time

# import the module from the build/ folder
sys.path.append('../build/')
from HHSP import HHSP



def main(dimension, dataset_size, radii, num_queries):

    # iniitialize the dataset: contiguous array, 4-byte floats
    data = np.ascontiguousarray(np.random.random([dataset_size,dimension]),dtype=np.float32)*2 - 1 # from [-1,1]
    print(data.flags)
    data_buffer = data.data

    # create testset
    query_data = np.ascontiguousarray(np.random.random([num_queries,dimension]),dtype=np.float32)*2 - 1 # from [-1,1]
    query_buffer = query_data.data

    # initialize the HHSP with dataset
    alg = HHSP(dimension)
    alg.add_dataset(data_buffer)

    # perform normal HSP Search
    start_time = time.time()
    gt_neighbors = alg.HSP_Test(query_buffer)
    gt_time = (time.time()-start_time)/num_queries
    print(f"Normal HSP Time (s): {gt_time:.4}")

    # create index
    alg.create_index(radii)
    alg.validate_index()

    # perform hierarchical HSP Search
    start_time = time.time()
    hhsp_neighbors = alg.Hierarchical_HSP_Test(query_buffer)
    hhsp_time = (time.time()-start_time)/num_queries
    print(f"Hierarchical HSP Time (s): {hhsp_time:.4}")
    print(f"Speedup over GT: {gt_time/hhsp_time:.4}")

    # compare the results
    for i in range(num_queries):
        if (gt_neighbors[i] != hhsp_neighbors[i]):
            print("Incorrect:")
            print("  - GT: ",gt_neighbors[i])
            print("  - HHSP: ",hhsp_neighbors[i])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", "-D",
        type=int,
        default=2
    )
    parser.add_argument("--size","-N",
        type=int,
        default=10000,
    )
    parser.add_argument('-r','--radii', 
        nargs='+',
        type=float,
        required=True)
    parser.add_argument(
        "--queries", "-T",
        type=int,
        default=100,
    )

    args = parser.parse_args()
    main(args.dimension, args.size, args.radii, args.queries)
