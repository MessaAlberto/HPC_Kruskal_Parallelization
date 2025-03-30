#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

typedef struct Edge {
  int u;
  int v;
  int weight;
} Edge;

typedef struct Subset {
  int parent;
  int rank;
} Subset;

void get_argv(int argc, char* argv[], int* V) {
  if (argc != 2) {
    printf("Usage: ./graph_generator <number_of_vertices>\n");
    exit(1);
  }

  *V = atoi(argv[1]);

  if (*V <= 1) {
    printf("Number of vertices must be greater than 1.\n");
    exit(1);
  }
}

void complete_graph(Edge** graph, int V, int max_wgt) {
  *graph = (Edge*)malloc(sizeof(Edge) * (int64_t)V * (V - 1) / 2);
  if (*graph == NULL) {
    printf("Memory allocation failed for complete graph.\n");
    exit(EXIT_FAILURE);
  }

  int64_t edge_count = 0;
  for (int i = 0; i < V; i++) {
    for (int j = i + 1; j < V; j++) {
      (*graph)[edge_count].u = i;
      (*graph)[edge_count].v = j;
      (*graph)[edge_count].weight = 1 + rand() % max_wgt;
      edge_count++;
    }
  }
}

// Calculate the number of missing edges to reach the desired density
int64_t calculate_missing_edges(int V, int density) {
  int64_t max_edges = (int64_t)V * (V - 1) / 2;
  float density_factor = density / 10.0;
  int64_t density_edges = (int64_t)(max_edges * density_factor);
  return density_edges - (int64_t)(V - 1);
}

void sparse_graph(Edge** graph, int V, int density, int max_wgt) {
  int64_t more_edges = calculate_missing_edges(V, density);

  *graph = (Edge*)malloc(sizeof(Edge) * (V - 1 + more_edges));
  int* nodes = (int*)malloc(sizeof(int) * V);
  uint8_t* adj_matrix = (uint8_t*)calloc((int64_t)V * V, sizeof(uint8_t));

  if (!(*graph) || !nodes || !adj_matrix) {
    printf("Memory allocation failed for sparse graph.\n");
    exit(EXIT_FAILURE);
  }

  // Initialize the nodes
  for (int i = 0; i < V; i++) {
    nodes[i] = i;
  }

  // Shuffle the nodes
  for (int i = 0; i < V; i++) {
    int j = rand() % V;
    int temp = nodes[i];
    nodes[i] = nodes[j];
    nodes[j] = temp;
  }
  
  // Connect all nodes
  int64_t edge_count = 0;
  for (int i = 1; i < V; i++) {
    int u = nodes[i];
    int v = nodes[rand() % i];
    int w = 1 + rand() % max_wgt;

    (*graph)[edge_count].u = u;
    (*graph)[edge_count].v = v;
    (*graph)[edge_count].weight = w;
    edge_count++;

    adj_matrix[(int64_t)u * V + v] = 1;
    adj_matrix[(int64_t)v * V + u] = 1;
  }  // All nodes are now connected

  // Add more edges to reach the desired density
  int64_t i = 0;
  while (i < more_edges) {
    int u = rand() % V;
    int v = rand() % V;
    int w = 1 + rand() % max_wgt;

    if (u != v && adj_matrix[(int64_t)u * V + v] == 0) {
      (*graph)[edge_count].u = u;
      (*graph)[edge_count].v = v;
      (*graph)[edge_count].weight = w;
      edge_count++;

      adj_matrix[(int64_t)u * V + v] = 1;
      adj_matrix[(int64_t)v * V + u] = 1;
      i++;
    }
  }

  free(adj_matrix);
  free(nodes);
}

void print_field(double times[10][10], int algo_times, const char* field_name) {
  printf("%s\n", field_name);
  printf(" nproc  ");
  for (int i = 0; i < 10; i++) {
    printf("| %6d%%    ", (i + 1) * 10);
  }
  printf("\n");

  for (int i = 0; i < algo_times; i++) {
    printf("  %3d   ", (int)pow(2, i));
    for (int j = 0; j < 10; j++) {
      printf("| %10.8f ", times[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_evaluation(double times[10][10], int algo_times) {
  printf("\n========= Speedup =========\n");
  printf(" thread  ");
  for (int i = 0; i < 10; i++) {
    printf("| %6d%%    ", (i + 1) * 10);
  }
  printf("\n");

  for (int i = 0; i < algo_times; i++) {
    printf("  %3d   ", (int)pow(2, i));
    for (int j = 0; j < 10; j++) {
      double speedup = times[0][j] / times[i][j];
      printf("| %10.8f ", speedup);
    }
    printf("\n");
  }
  printf("\n");

  printf("\n========= Efficiency =========\n");
  printf(" thread  ");
  for (int i = 0; i < 10; i++) {
    printf("| %6d%%    ", (i + 1) * 10);
  }
  printf("\n");

  for (int i = 0; i < algo_times; i++) {
    printf("  %3d   ", (int)pow(2, i));
    for (int j = 0; j < 10; j++) {
      double efficiency = (times[0][j] / times[i][j]) / pow(2, i);
      printf("| %10.8f ", efficiency);
    }
    printf("\n");
  }
  printf("\n");
}

int cmpfunc(const void* a, const void* b) {
  const Edge* edgeA = (const Edge*)a;
  const Edge* edgeB = (const Edge*)b;
  if (edgeA->weight < edgeB->weight) return -1;
  if (edgeA->weight > edgeB->weight) return 1;
  return 0;
}

void merge(Edge* graph, int64_t l, int64_t m, int64_t r) {
  Edge* tmp = (Edge*)malloc(sizeof(Edge) * (r - l + 1));
  if (!tmp) {
    printf("Memory allocation failed for merge.\n");
    exit(EXIT_FAILURE);
  }

  int64_t i = l, j = m + 1, k = 0;

  while (i <= m && j <= r) {
    tmp[k++] = (graph[i].weight <= graph[j].weight) ? graph[i++] : graph[j++];
  }

  while (i <= m) tmp[k++] = graph[i++];
  while (j <= r) tmp[k++] = graph[j++];

  for (i = l, k = 0; i <= r; i++, k++) {
    graph[i] = tmp[k];
  }

  free(tmp);
}

void mergeSortParallel(Edge* graph, int64_t E) {
  int64_t curr_size, start;
  int num_threads = omp_get_max_threads();
  int64_t block_size = E / num_threads;

  #pragma omp parallel for
  for (int64_t i = 0; i < num_threads; i++) {
    int64_t start = i * block_size;
    int64_t end = (i == num_threads - 1) ? E - 1 : start + block_size - 1;
    qsort(graph + start, (size_t)(end - start + 1), sizeof(Edge), cmpfunc);
  }

  for (curr_size = block_size; curr_size <= E; curr_size *= 2) {
    #pragma omp parallel for private(start)
    for (start = 0; start < E; start += 2 * curr_size) {
      int64_t mid = start + curr_size - 1;
      int64_t end = (start + 2 * curr_size - 1 < E) ? start + 2 * curr_size - 1 : E - 1;

      if (mid < E && graph[mid].weight > graph[mid + 1].weight) {
        merge(graph, start, mid, end);
      }
    }
  }
}

int main(int argc, char* argv[]) {

  int nproc = omp_get_max_threads();
  printf("Max number of threads: %d\n", nproc);

  int V = 0;
  int64_t E = 0;
  Edge* graph = NULL;
  Edge* original_graph = NULL;

  // Calculate the number of times the algorithm will run
  // Sequential Kruskal + Parallel Kruskal with 2, 4, 8, ... processes
  int algo_times = 1 + (int)log2(nproc);
  double times[algo_times][10];   // 10 graph densities from 10% to 100%

  double start, end;

  get_argv(argc, argv, &V);
  int max_wgt = INT32_MAX / (V - 1);

  srand(time(NULL));

  // Graph density from 10% to 100%
  for (int i = 1; i <= 10; i++) {    
    printf("===== Density: %d%% =====\n", i * 10);
    fflush(stdout);

    if (i < 10) {
      E = V - 1 + calculate_missing_edges(V, i);
      if (E <= V - 1)
        continue;  // Skip densities with less than V - 1 edges (not connected graph)
      sparse_graph(&graph, V, i, max_wgt);
    } else {
      E = (int64_t)V * (V - 1) / 2;
      complete_graph(&graph, V, max_wgt);
    }

    // Save the original graph
    original_graph = (Edge*)malloc(sizeof(Edge) * E);
    if (original_graph == NULL) {
      printf("Memory allocation failed for original graph.\n");
      return 1;
    }

    memcpy(original_graph, graph, sizeof(Edge) * E);
    
    // Iterate the number of processes, increasing by powers of 2 (1, 2, 4, 8, ...)
    for (int j = 1; j <= nproc; j *= 2) {
      omp_set_num_threads(j);
      
      printf("Parallel Kruskal: %d threads ...\n", j);
      fflush(stdout);

      // Restore the original graph
      memcpy(graph, original_graph, sizeof(Edge) * E);

      printf("Let's do sorting...\n");
      fflush(stdout);
    
      start = omp_get_wtime();
      mergeSortParallel(graph, E);
      end = omp_get_wtime();

      for (int64_t k = 1; k < E; k++) {
        if (graph[k].weight < graph[k - 1].weight) {
          printf("Error: Parallel sorting failed.\n");
          return 1;
        }
      }

      printf("Graph sorted\n");
      fflush(stdout);

      int run_index = (int)log2(j);
      times[run_index][i - 1] = end - start;

      printf("Done\n\n");
      fflush(stdout);
    }

    free(original_graph);
    free(graph);
  }

  // Print the results
  printf("\n\n");
  printf(" ========== Results ==========\n");

  print_field(times, algo_times, "Sorting Time:");
  print_evaluation(times, algo_times);

  return 0;
}