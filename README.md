# HPC Kruskal Parallelization

## How to Run a Test

### 1. Graph Generation

#### Go under the utils folder
```bash
cd utils/
```

#### Compile files
```bash
make
```

#### Create a connected graph with a random number of edges
```bash
./graph_gen <number_of_nodes> <max_weight>
```
- This command will generate a graph with a random number of edges and weights between 1 and max_weight.

Example: 
```bash
./graph_gen 10000 200
```
Output:
```
Graph with 10000 vertices and 0.96 density generated.
```
File `graph_10000_nodes_0.96_density.txt` will be created in the folder `graphs/`.

#### Create a complete graph
```bash
./complete_graph <number_of_nodes> <max_weight>
```

Example: 
```bash
./complete_graph 10000 200
```
Output:
```
Complete graph with 10000 vertices generated.

```
File `graph_10000_nodes_1.00_density.txt` will be created in the folder `graphs/`.


### 2. Kruskal Algorithm
Go to src folder
```bash
cd src/
```

Compile files
```bash
make
```

#### Run the sequential version
```bash
./seq_kruskal <input_file>
```
Example:
```bash
./seq_kruskal ../graphs/graph_10000_nodes_0.96_density.txt
```
Output saved in `../mst/seq_mst_graph_10000_nodes_0.96_density.txt`.

#### Run the parallel version
```bash
./mpiexec -n <number_of_processes> ./par_kruskal <input_file>
```
Example:
```bash
./mpiexec -n 4 ./par_kruskal ../graphs/graph_10000_nodes_0.96_density.txt
```
Output saved in `../mst/mpi_mst_graph_10000_nodes_0.96_density.txt`.

### 3. Check if the MST is correct
(Still to implement automatic check)
#### Manual check
Open files:
- `../mst/seq_mst_graph_10000_nodes_0.96_density.txt`
- `../mst/mpi_mst_graph_10000_nodes_0.96_density.txt`
First line of the file has the format: `<number_of_nodes> <number_of_edges> <total_weight>`.
Next lines have the format: `<node1> <node2> <weight>`.

Check if mst are connected:
```bash
cd ../utils/
```
```bash
./check_connectivity <input_file>
```
Example:
```bash
./check_connectivity ../mst/seq_mst_graph_10000_nodes_0.96_density.txt
```
and
```bash
./check_connectivity ../mst/mpi_mst_graph_10000_nodes_0.96_density.txt
```




