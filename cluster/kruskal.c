#include <inttypes.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct Edge {
  int u;
  int v;
  int weight;
} Edge;

typedef struct Subset {
  int parent;
  int rank;
} Subset;

typedef struct AlgoTimes {
  double qsort_time;
  double distribute_time;
  double kruskal_time;
  double collect_time;
  double total_time;
} AlgoTimes;

typedef double (*FieldGetter)(AlgoTimes*);

double get_qsort_time(AlgoTimes* times) { return times->qsort_time; }

double get_distribute_time(AlgoTimes* times) { return times->distribute_time; }

double get_kruskal_time(AlgoTimes* times) { return times->kruskal_time; }

double get_collect_time(AlgoTimes* times) { return times->collect_time; }

double get_total_time(AlgoTimes* times) { return times->total_time; }

MPI_Datatype MPI_EDGE;

int compare_edges(const void* a, const void* b) {
  return ((Edge*)a)->weight - ((Edge*)b)->weight;
}

// Find the root of the set
int find(Subset* sets, int i) {
  if (sets[i].parent != i) {
    sets[i].parent = find(sets, sets[i].parent);
  }
  return sets[i].parent;
}

// Union two sets
void union_sets(Subset* sets, int root_u, int root_v) {
  if (sets[root_u].rank < sets[root_v].rank) {
    sets[root_u].parent = root_v;
  } else if (sets[root_u].rank > sets[root_v].rank) {
    sets[root_v].parent = root_u;
  } else {
    sets[root_v].parent = root_u;
    sets[root_u].rank++;
  }
}

void create_mpi_edge_type() {
  int block_lengths[3] = {1, 1, 1};
  MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint offsets[3];

  offsets[0] = offsetof(Edge, u);
  offsets[1] = offsetof(Edge, v);
  offsets[2] = offsetof(Edge, weight);

  MPI_Type_create_struct(3, block_lengths, offsets, types, &MPI_EDGE);
  MPI_Type_commit(&MPI_EDGE);
}

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
  *graph = (Edge*)malloc(sizeof(Edge) * V * (V - 1) / 2);
  if (*graph == NULL) {
    printf("Memory allocation failed for complete graph.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
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
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 1);
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

    if (u != v && adj_matrix[u * V + v] == 0) {
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

void kruskal(Edge* graph, int V, int64_t E, Edge** mst, int* mst_edges, int* mst_weight,
             Subset* sets) {
  for (int i = 0; i < V; i++) {
    sets[i].parent = i;
    sets[i].rank = 0;
  }

  // MST data structure
  *mst = (Edge*)malloc(sizeof(Edge) * (V - 1));
  if (*mst == NULL) {
    printf("Memory allocation failed for MST.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  *mst_weight = 0;
  *mst_edges = 0;

  for (int64_t i = 0; i < E; i++) {
    int u = graph[i].u;
    int v = graph[i].v;

    int root_u = find(sets, u);
    int root_v = find(sets, v);

    if (root_u != root_v) {
      // Include the edge in the MST
      (*mst)[(*mst_edges)++] = graph[i];
      *mst_weight += graph[i].weight;

      // Union the sets
      union_sets(sets, root_u, root_v);
    }
  }
}

void distibute_data(Edge* graph, Edge** local_graph, int64_t E, int* send_counts, int* displs,
                    int nproc, int rank, MPI_Comm active_comm) {
  if (E < nproc) {
    printf("Error: more processes than edges.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int64_t base_size = E / nproc;
  int64_t extra = E % nproc;
  if (base_size + extra > INT_MAX) {
    printf("Error: data is too large to distribute.\n");
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int offset = 0;
  for (int i = 0; i < nproc; i++) {
    send_counts[i] = base_size + (i < extra ? 1 : 0);
    displs[i] = offset;
    offset += send_counts[i];
  }

  int local_E = send_counts[rank];

  *local_graph = (Edge*)malloc(sizeof(Edge) * local_E);
  if (*local_graph == NULL) {
    printf("Memory allocation failed for local graph.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Scatterv(graph, send_counts, displs, MPI_EDGE, *local_graph, local_E, MPI_EDGE, 0,
               active_comm);

  if (*local_graph == NULL) {
    printf("Error distributing data.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void collect_mst(int rank, int nproc, int V, Edge** local_mst, int* mst_edges, int* mst_weight,
                 Subset* sets, MPI_Comm active_comm) {
  int power = 1;
  while (nproc > 1) {
    if ((rank / power) % 2 == 0) {
      int recv_rank = rank + power;
      int recv_edges = 0;
      MPI_Recv(&recv_edges, 1, MPI_INT, recv_rank, 0, active_comm, MPI_STATUS_IGNORE);

      Edge* recv_mst = (Edge*)calloc(recv_edges, sizeof(Edge));

      if (recv_mst == NULL) {
        printf("Memory allocation failed for received MST.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      MPI_Recv(recv_mst, recv_edges, MPI_EDGE, recv_rank, 0, active_comm, MPI_STATUS_IGNORE);

      for (int j = 0; *mst_edges < V - 1 && j < recv_edges; j++) {
        int u = recv_mst[j].u;
        int v = recv_mst[j].v;

        int root_u = find(sets, u);
        int root_v = find(sets, v);

        if (root_u != root_v) {
          (*local_mst)[(*mst_edges)++] = recv_mst[j];
          *mst_weight += recv_mst[j].weight;
          union_sets(sets, root_u, root_v);
        }
      }

      free(recv_mst);
    } else {
      int send_rank = rank - power;
      MPI_Send(mst_edges, 1, MPI_INT, send_rank, 0, active_comm);
      MPI_Send(*local_mst, *mst_edges, MPI_EDGE, send_rank, 0, active_comm);

      break;
    }
    power *= 2;
    nproc /= 2;
  }
}

void print_field(AlgoTimes times[10][10], int algo_times, const char* field_name,
                 FieldGetter field) {
  printf("%s\n", field_name);
  printf(" nproc  ");
  for (int i = 0; i < 10; i++) {
    printf("| %6d%%    ", (i + 1) * 10);
  }
  printf("\n");

  for (int i = 0; i < algo_times; i++) {
    printf("  %3d   ", (int)pow(2, i));
    for (int j = 0; j < 10; j++) {
      printf("| %10.8f ", field(&times[i][j]));
    }
    printf("\n");
  }
  printf("\n");
}

void print_evaluation(AlgoTimes times[10][10], int algo_times) {
  printf("\n========= Task Evaluation on complete graph (%% of Total Time) =========\n");
  printf(" nproc\n");

  for (int i = 0; i < algo_times; i++) {
    printf("  %3d   ", (int)pow(2, i));
    double total = times[i][9].total_time;

    double qsort_percentage = (times[i][9].qsort_time / total) * 100.0;
    double distribute_percentage = (times[i][9].distribute_time / total) * 100.0;
    double kruskal_percentage = (times[i][9].kruskal_time / total) * 100.0;
    double collect_percentage = (times[i][9].collect_time / total) * 100.0;

    printf("| Qsort: %5.1f%% \tDist: %5.1f%% \tKruskal: %5.1f%% \tCollect: %5.1f%%\n",
           qsort_percentage, distribute_percentage, kruskal_percentage, collect_percentage);
  }
  printf("\n========= Speedup =========\n");
  printf(" nproc  ");
  for (int i = 0; i < 10; i++) {
    printf("| %6d%%    ", (i + 1) * 10);
  }
  printf("\n");

  for (int i = 0; i < algo_times; i++) {
    printf("  %3d   ", (int)pow(2, i));
    for (int j = 0; j < 10; j++) {
      double speedup = times[0][j].total_time / times[i][j].total_time;
      printf("| %10.8f ", speedup);
    }
    printf("\n");
  }
  printf("\n");

  printf("\n========= Efficiency =========\n");
  printf(" nproc  ");
  for (int i = 0; i < 10; i++) {
    printf("| %6d%%    ", (i + 1) * 10);
  }
  printf("\n");

  for (int i = 0; i < algo_times; i++) {
    printf("  %3d   ", (int)pow(2, i));
    for (int j = 0; j < 10; j++) {
      double efficiency = (times[0][j].total_time / times[i][j].total_time) / pow(2, i);
      printf("| %10.8f ", efficiency);
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  create_mpi_edge_type();
  int V = 0;
  int64_t E = 0;
  int mst_weight = 0, mst_edges = 0;
  Edge* graph = NULL;
  Edge* original_graph = NULL;
  Edge* mst = NULL;
  Subset* sets = NULL;

  int send_counts[nproc], displs[nproc];
  Edge* local_graph = NULL;
  Edge* local_mst = NULL;

  int result_mst_wgt = 0;
  // Calculate the number of times the algorithm will run
  // Sequential Kruskal + Parallel Kruskal with 2, 4, 8, ... processes
  int algo_times = 1 + (int)log2(nproc);
  AlgoTimes times[algo_times][10];  // 10 graph densities from 10% to 100%

  double start_time, end_time;
  double qsort_start_time, qsort_end_time;
  double distribute_start_time, distribute_end_time;
  double kruskal_start_time, kruskal_end_time;
  double collect_start_time, collect_end_time;

  get_argv(argc, argv, &V);
  int max_wgt = INT32_MAX / (V - 1);

  srand(time(NULL));

  // Graph density from 10% to 100%
  for (int i = 1; i <= 10; i++) {
    if (rank == 0) {
      printf("----- Density: %d%% -----\n", i * 10);
      fflush(stdout);
    }

    if (i < 10) {
      E = V - 1 + calculate_missing_edges(V, i);
      if (E <= V - 1) continue;  // Skip densities with less than V - 1 edges

      if (rank == 0) {  // Process 0 generates the graph
        sparse_graph(&graph, V, i, max_wgt);
      }
    } else {
      E = V * (V - 1) / 2;
      if (rank == 0) {  // Process 0 generates the graph
        complete_graph(&graph, V, max_wgt);
      }
    }

    // Save the original graph
    if (rank == 0) {
      original_graph = (Edge*)malloc(sizeof(Edge) * E);
      if (original_graph == NULL) {
        printf("Memory allocation failed for original graph.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      memcpy(original_graph, graph, sizeof(Edge) * E);
    }

    // Broadcast the number of vertices
    MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // =========== Sequential Kruskal Algorithm ===========
    if (rank == 0) {
      printf("Sequential Kruskal: ... ");
      fflush(stdout);
      sets = (Subset*)malloc(sizeof(Subset) * V);
      if (sets == NULL) {
        printf("Memory allocation failed for sets.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      start_time = MPI_Wtime();
      qsort_start_time = start_time;

      qsort(graph, E, sizeof(Edge), compare_edges);

      qsort_end_time = MPI_Wtime();
      kruskal_start_time = qsort_end_time;

      kruskal(graph, V, E, &mst, &mst_edges, &mst_weight, sets);

      kruskal_end_time = MPI_Wtime();
      end_time = kruskal_end_time;

      result_mst_wgt = mst_weight;
      times[0][i - 1].qsort_time = qsort_end_time - qsort_start_time;
      times[0][i - 1].distribute_time = 0;
      times[0][i - 1].kruskal_time = kruskal_end_time - kruskal_start_time;
      times[0][i - 1].collect_time = 0;
      times[0][i - 1].total_time = end_time - start_time;

      printf("Done\n");
      fflush(stdout);
    }

    // =========== Parallel Kruskal Algorithm ===========
    // Iterate the number of processes, increasing by powers of 2 (2, 4, 8, ...)
    for (int j = 2; j <= nproc; j *= 2) {
      if (rank == 0) {
        printf("Parallel Kruskal: %d processes ... ", j);
        fflush(stdout);

        // Restore the original graph
        memcpy(graph, original_graph, sizeof(Edge) * E);
      }

      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Comm active_comm = MPI_COMM_NULL;
      MPI_Group world_group, active_group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);

      if (rank < j) {  // Create the active group
        int ranks[j];
        for (int k = 0; k < j; k++) ranks[k] = k;
        MPI_Group_incl(world_group, j, ranks, &active_group);
      } else {
        MPI_Group_excl(world_group, 1, &rank, &active_group);
      }

      MPI_Comm_create(MPI_COMM_WORLD, active_group, &active_comm);

      if (active_comm == MPI_COMM_NULL) {  // Processes not in the active group
        continue;
      }

      sets = (Subset*)malloc(sizeof(Subset) * V);
      if (sets == NULL) {
        printf("Memory allocation failed for sets.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      if (rank == 0) {
        start_time = MPI_Wtime();
        qsort_start_time = start_time;

        qsort(graph, E, sizeof(Edge), compare_edges);

        qsort_end_time = MPI_Wtime();
        distribute_start_time = qsort_end_time;
      }

      // Distribute the data to the processes
      distibute_data(graph, &local_graph, E, send_counts, displs, j, rank, active_comm);

      if (rank == 0) {
        distribute_end_time = MPI_Wtime();
        kruskal_start_time = distribute_end_time;
      }

      // Compute the local MST
      kruskal(local_graph, V, send_counts[rank], &local_mst, &mst_edges, &mst_weight, sets);

      if (rank == 0) {
        kruskal_end_time = MPI_Wtime();
        collect_start_time = kruskal_end_time;
      }

      // Collect the MST from all the processes
      collect_mst(rank, j, V, &local_mst, &mst_edges, &mst_weight, sets, active_comm);

      if (rank == 0) {
        collect_end_time = MPI_Wtime();
        end_time = collect_end_time;

        int run_index = (int)log2(j);
        times[run_index][i - 1].qsort_time = qsort_end_time - qsort_start_time;
        times[run_index][i - 1].distribute_time = distribute_end_time - distribute_start_time;
        times[run_index][i - 1].kruskal_time = kruskal_end_time - kruskal_start_time;
        times[run_index][i - 1].collect_time = collect_end_time - collect_start_time;
        times[run_index][i - 1].total_time = end_time - start_time;

        if (mst_weight != result_mst_wgt) {
          printf("Error: MST weight is different from sequential.\n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        printf("Done\n");
        fflush(stdout);
      }
    }
  }

  // Print the results
  if (rank == 0) {
    printf("\n\n");
    printf(" ========== Results ==========\n");

    print_field(times, algo_times, "Qsort Time:", get_qsort_time);
    print_field(times, algo_times, "Distribute Time:", get_distribute_time);
    print_field(times, algo_times, "Kruskal Time:", get_kruskal_time);
    print_field(times, algo_times, "Collect Time:", get_collect_time);
    print_field(times, algo_times, "Total Time:", get_total_time);

    print_evaluation(times, algo_times);
  }

  free(graph);
  free(original_graph);
  free(mst);
  free(sets);
  free(local_graph);
  free(local_mst);

  MPI_Type_free(&MPI_EDGE);
  MPI_Finalize();

  return 0;
}