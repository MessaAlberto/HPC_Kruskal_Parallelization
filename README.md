# HPC Kruskal Parallelization Project

This repository contains the implementation of a parallelized version of **Kruskal's algorithm** for finding the Minimum Spanning Tree (MST) of a graph. The project is designed to run on both local and cluster environments, utilizing MPI for parallelization.

## Sequential Kruskal Algorithm
### Overview
Kruskal's algorithm is a **greedy algorithm** used to find the **Minimum Spanning Tree (MST)** of a connected, undirected, and weighted graph. It works by sorting all the edges in increasing order of their weights and adding them one by one to the MST, **as long as they don’t form a cycle**.

To efficiently check whether an edge forms a cycle, it uses the **Disjoint Set Union (DSU)** data structure (also known as Union-Find).

---

### Steps

1. Order the edges by increasing weight.
2. Initialize a DSU structure with all vertices in separate sets.
3. Iterate through the sorted edges and:
   - If the current edge connects two different components, add it to the MST and unite the components.
   - Otherwise, skip it (would form a cycle).
4. Repeat until the MST contains exactly `V - 1` edges (where `V` is the number of vertices).

---

### Pseudocode

```c
// graph: array of edges
// V: number of vertices
// E: number of edges

sort(graph);  // O(E log E)
init_dsu(V);  // O(V)

MST = []

for (i = 0; i < E && MST.size() < V - 1; i++) {
    u = graph[i].u
    v = graph[i].v

    if (find(u) != find(v)) {
        MST.push(graph[i]);  // add edge to MST
        union_sets(u, v);  // unite sets
    }
}
```

## Parallel Kruskal Algorithm – Idea and Approach

This project explores the parallelization of Kruskal’s algorithm by focusing on the most computationally expensive phase: **edge sorting**. We implemented two distinct MPI-based approaches to parallel sorting, each applied within the context of Kruskal's algorithm for computing the Minimum Spanning Tree (MST) of a graph.

### Parallel Strategies

The core idea behind our parallel approach is to **distribute edges by weight** across multiple processes. Edges with the lowest weights are assigned to the lowest-ranked process, the next lightest edges to the next process, and so forth. This ensures that each process holds a subset of edges with lighter weights than those held by higher-ranked processes.

Each process then independently computes a **local MST** from its assigned edges. Importantly, the MST built by a lower-ranked process will always contain lighter edges than that of any higher-ranked process. During the merging phase, a process receives the MST from the next higher-ranked process and adds only the missing edges to its own MST—**no additional sorting is needed**.

- **No re-sorting** is required after merging.
- **Global edge order is preserved** through careful distribution.
- **Efficient merging** is achieved by avoiding redundant comparisons and edge duplication.

This design avoids repeated sorting of the same edges, at the cost of a controlled overhead for distributing edges. We provide two approaches to perform this distribution efficiently:

### 1. `pivotSort` – Parallel Quicksort Approach
Assuming a near-uniform distribution of edge weights (supported by the central limit theorem for large datasets), we divide the global weight range into equal intervals. This allows us to assign edges to processes using a pivoting strategy:
1. **Pivot computation**: Process 0 determines pivot values by dividing the maximum edge weight range into intervals based on the number of processes.
2. **Edge assignment**: Process 0 iterates through all edges, assigning each to a process based on which pivot interval it falls into.
3. **Data distribution**: Each process receives a chunk of edges corresponding to its weight interval.

After receiving their chunk, each process sorts its local edges and runs Kruskal's algorithm independently on its local graph.

### 2. `mergeSort` – Global Parallel Mergesort
The second strategy is based on a parallel mergesort-like approach. The key idea is to **first sort the entire edge set in parallel**, then redistribute the sorted edges:
 1. **Initial distribution**: Process 0 splits the global edge list and sends chunks to all processes.
 2. **Local sorting**: Each process sorts its own subset of edges independently (e.g., using quicksort).
 3. **Parallel merging**: Processes merge their sorted chunks in a `tree-like hierarchy` to minimize communication steps. For example:
    ```bash
     rank0  rank1  rank2  rank3
         \  /          \  /
          \/            \/
          rank0      rank1
               \    /
                \  /
                rank0
    ```
 4. **Global collection**: Process 0 eventually gathers the fully sorted edge list.
 5. **Final redistribution**: The sorted edges are split and reassigned to processes based on rank, with lighter edges going to lower ranks.

Once again, each process runs Kruskal’s algorithm locally on its assigned edges.

### Trade-off between `pivotSort` and `mergeSort`-like Method

In the **`pivotSort`** approach, **Process 0** is responsible for iterating through the entire graph, filling each process’s chunk with edges based on the computed pivot intervals. Afterward, it sends the respective chunk to each process. While this approach ensures that each process gets a chunk of edges, it requires **one round of communication** and **less overhead** on the initial data distribution.

On the other hand, the **`mergeSort`-like approach** distributes the graph multiple times:
1. Initially, **Process 0** distributes the unordered graph to all processes.
2. Each process sorts its edges locally.
3. Afterward, a **tree-like merge** occurs to combine the sorted chunks, which introduces **additional communication rounds**.
4. Finally, **Process 0** collects the fully sorted graph and redistributes the edges, assigning chunks based on process ranks, ensuring that the lower-ranked processes receive lighter edges.

Interestingly, despite the apparent differences in distribution and communication complexity, **both approaches yield similar runtimes**:

- **`pivotSort`** is faster due to its lower communication overhead, while
- **`mergeSort`** offers better load balancing but incurs higher communication costs.

In practice, however, both approaches converge to **comparable runtimes**, as they strike a balance between computation and communication.

---

### MST Merging in Parallel Kruskal

In both sorting algorithms, each process ends up with a **local MST** that contains edges with smaller weights compared to the MSTs of higher-ranked processes. This is crucial because it enables an **efficient merge process**. When a process receives an MST from a higher-ranked process, it knows that its MST is lighter and thus does not need to be checked again. Instead, the process simply adds missing edges from the incoming MST.

This eliminates the need for additional sorting or adjustments, as the lower-ranked process is guaranteed to have the lighter MST. This results in a **faster merge** without redundant operations, ensuring that the final MST is built efficiently.

## Project Structure
```
├── cluster
│   ├── Makefile
│   ├── output/
│   ├── output_err/
│   ├── kruskal_mergeSort.c
│   ├── kruskal_pivotSort.c
│   └── script/
│       ├── qsub_dir.sh
│       └── test.sh
└── local
    ├── graphs/
    ├── mst/
    ├── kruskal_algorithm/
    │   ├── Makefile
    |   ├── include/
    │   │   ├── disjoint_set.h
    │   │   └── graph_utils.h
    │   ├── src/
    │   │   ├── disjoint_set.c
    │   │   └── graph_utils.c
    │   ├── par_kruskal.c
    │   └── seq_kruskal.c
    └── utils/
        ├── Makefile
        ├── include/graph_gen_utils.h
        ├── src/graph_gen_utils.c
        ├── check_connectivity.c
        ├── complete_graph_gen.c
        └── sparse_graph_gen.c
```


## Cluster
The `cluster/` folder contains the necessary files to run the Kruskal algorithm on a cluster using MPI. A single job executes the Kruskal algorithm with the code in `kruskal_mergeSort.c` or `kruskal_pivotSort.c` on a graph of size `V`, where `V` represents the number of nodes in the graph. The program is designed to run with the same number of nodes and process multiple graphs of increasing size and density.

In `kruskal_mergeSort.c` and `kruskal_pivotSort.c`, there are two `#define` statements:
```c
#define DENSITY_START 1  // 10% density
#define DENSITY_END 10   // 100% density
```
This configuration means that the program will run on graphs with densities ranging from 10% to 100%, increasing by 10% at each step. Therefore, the program will execute on graphs of size V with 10%, 20%, ..., up to 100% density.

Additionally, depending on the number of MPI processes, the program will run the algorithm using 2, 4, 8, ..., up to the total number of available processes.

This allows us to gain a comprehensive overview of the results across different graph sizes and varying numbers of MPI processes.

### How to Run a Test
1. **Select the algorithm:**

   Choose between `kruskal_mergeSort.c` or `kruskal_pivotSort.c` by editing the `Makefile`. Uncomment one of the following lines:

     ```makefile
     MAIN_SRC = kruskal_pivotSort.c
     # MAIN_SRC = kruskal_mergeSort.c
     ```

2. **Set the test parameters:**
   Navigate to the `script/` folder and edit the `test.sh` file. The script will submit multiple PBS jobs to excecute the Kruskal algorithm on graphs of increasing size by power of 2. Change the following parameters to submit jobs for different graphs:

   ```bash
   # Test graph of 1024 to 65536 nodes
   v_start=1024
   v_end=65536
   ```
   This will run the algorithm on graphs with 1024, 2048, 4096, ..., up to 65536 nodes.

   You can custom output file names. Uncomment the desired line and adjust the output paths as needed:
   ```bash
   qsub -v V="$v" -o "./output/pivotSort_${v}_nodes.out" -e "./output_err/pivotSort_${v}_nodes.err" $SCRIPT_NAME
   # qsub -v V="$v" -o "./output/mergeSort_${v}_nodes.out" -e "./output_err/mergeSort_${v}_nodes.err" $SCRIPT_NAME
   ```

3. **Configure job parameters:**
   Edit the `qsub_dir.sh` file to set the number of nodes, processors and memory requirements per node. The default values are set to 8 nodes, 64 processors per node, and 128gb of memory in order to run the kruskal parallel version with 512 processes. Adjust these values according to your cluster configuration.

   ```bash
    #PBS -l select=8:ncpus=64:mem=128gb
   ```
   If you're increasing graph size significantly, it's recommended to increase the memory to `512gb` due to the worst-case memory usage by `rank 0`, which may need to store up to three full graph copies.

4. **Run the test:**
   To execute the test, navigate to the root `cluster/` directory and use the `make` command. This will compile the source and launch the PBS jobs.

   ```bash
    cd cluster/
    make
   ```
   The job outputs and error logs will be saved in the `output/` and `output_err/` folders respectively.

## Local
The `local/` folder contains the tools needed to generate graphs and test the Kruskal algorithm, both in sequential and parallel versions. The steps to run the tests in a local environment are outlined below.

### How to Run a Test
1. **Navigate to the `utils` folder:**

   Go to the `local/utils` folder, where the programs to generate graphs and check connectivity are located.

   ```bash
   cd local/utils/
   ```

2. **Compile the files:**

   Run the `make` command to compile the necessary programs. Executables will be generated in the `bin` folder.

   ```bash
   make
   ```

   This command compiles the source files and creates three executables:
   - `sparse_graph_gen`
   - `complete_graph_gen`
   - `check_connectivity`

3. **Generate a complete graph:**

   Use the `complete_graph` executable to generate a complete graph with a specific number of vertices (e.g., 2000).

   ```bash
   ./bin/complete_graph_gen.exe 2000
   ```

4. **Navigate to the `kruskal_algorithm` folder:**

   Once the graph is generated, move to the `local/kruskal_algorithm` folder, which contains the Kruskal algorithm code.

   ```bash
   cd ../kruskal_algorithm/
   ```
5. **Makefile Configuration**

   Before compiling, **update the MPI paths** in the Makefile to match your system configuration:

   ```makefile
   MPIFLAGS := -I "your_path/MPI/Include"
   MPILIB := -L "your_path/MPI/Lib/x64" -lmsmpi
   ```
   Run `make` to compile both the sequential and parallel versions of the Kruskal algorithm.

   ```bash
   make
   ```

6. **Run the test:**

   Run tests on the Kruskal algorithm using the `make test` command, passing the path to the graph as a parameter. For example:

   ```bash
   make test GRAPH=../graphs/graph_2000_nodes_1.00_density.txt
   ```

   This command runs the sequential version followed by the parallel version (with MPI), showing execution times.

---