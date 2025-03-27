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
  double distribute_unordered_graph_time;
  double qsort_time;
  double merge_graph_time;
  double distribute_ordered_graph_time;
  double kruskal_time;
  double collect_mst_time;
  double total_time;
} AlgoTimes;

typedef double (*FieldGetter)(AlgoTimes*);

double get_distribute_unordered_graph_time(AlgoTimes* times) {
  return times->distribute_unordered_graph_time;
}

double get_qsort_time(AlgoTimes* times) { return times->qsort_time; }

double get_merge_graph_time(AlgoTimes* times) { return times->merge_graph_time; }

double get_distribute_ordered_graph_time(AlgoTimes* times) {
  return times->distribute_ordered_graph_time;
}

double get_kruskal_time(AlgoTimes* times) { return times->kruskal_time; }

double get_collect_mst_time(AlgoTimes* times) { return times->collect_mst_time; }

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
  *graph = (Edge*)malloc(sizeof(Edge) * (int64_t)V * (V - 1) / 2);
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

void distibute_data(Edge* graph, Edge** local_graph, int64_t E, int64_t* send_counts,
                    int64_t* displs, int nproc, int rank, MPI_Comm active_comm) {
  if (E < nproc) {
    printf("Error: more processes than edges.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int64_t base_size = E / nproc;
  int64_t extra = E % nproc;
  int64_t offset = 0;

  for (int i = 0; i < nproc; i++) {
    send_counts[i] = base_size + (int64_t)(i < extra ? 1 : 0);
    displs[i] = offset;
    offset += send_counts[i];
  }

  *local_graph = (Edge*)malloc(sizeof(Edge) * send_counts[rank]);
  // Send counts and displacements must be int
  int* send_counts_int = (int*)malloc(sizeof(int) * nproc);
  int* displs_int = (int*)malloc(sizeof(int) * nproc);

  if (!(*local_graph) || !send_counts_int || !displs_int) {
    printf("Memory allocation failed for local graph.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (E > INT_MAX) {
    if (rank == 0) {
      printf("Distribution in multiple batches...\n");
      fflush(stdout);
    }
  
    // Max number of edges to send to each process for each MPI_Scatterv
    int batch_size = INT_MAX / nproc;
    // Number of batches
    int64_t num_batches =
        ((int64_t)E + batch_size * nproc - 1) / (int64_t)(batch_size * nproc);
    // Max number of edges to send to each process in total
    int64_t offset = ((int64_t)E + nproc - 1) / (int64_t)nproc;

    Edge* ordered_graph = NULL;
    int64_t* batch_index = NULL;

    // Only rank 0 creates the ordered graph and batch index to distribute the data
    if (rank == 0) {
      ordered_graph = (Edge*)malloc(sizeof(Edge) * E);
      batch_index = (int64_t*)malloc(sizeof(int64_t) * nproc);

      if (!ordered_graph || !batch_index) {
        printf("Memory allocation failed for ordered graph.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      for (int i = 0; i < nproc; i++) {
        batch_index[i] = (int64_t)base_size + (i < extra ? 1 : 0) + i * offset;
      }

      int64_t count = 0;
      for (int batch = 0; batch < num_batches; batch++) {
        for (int proc = 0; proc < nproc; proc++) {
          int64_t start_pos = (int64_t)batch * batch_size + proc * offset;
          int64_t end_pos = start_pos + batch_size;

          if (end_pos > batch_index[proc]) {
            end_pos = batch_index[proc];
          }

          for (int64_t pos = start_pos; pos < end_pos && pos < E; pos++) {
            ordered_graph[count++] = graph[pos];
          }
        }
      }
    }

    // Due to MPI_Scatterv limitations(send_count must be int), we need to send the data in multiple batches
    for (int batch = 0; batch < num_batches; batch++) {
      int64_t start_pos = (int64_t)batch * batch_size * nproc;
      int64_t end_pos = start_pos + batch_size * nproc;
      if (end_pos > E) {
        end_pos = E;
      }

      int chunk_size = (int)(end_pos - start_pos);

      for (int i = 0; i < nproc; i++) {
        send_counts_int[i] = chunk_size / nproc + (i < chunk_size % nproc ? 1 : 0);
        displs_int[i] = i * send_counts_int[i];
      }

      MPI_Scatterv(ordered_graph + start_pos, send_counts_int, displs_int, MPI_EDGE,
                   *local_graph + batch_size * (int64_t)batch, send_counts_int[rank], MPI_EDGE, 0, active_comm);

      if (*local_graph == NULL) {
        printf("Error distributing data.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    free(ordered_graph);
    free(batch_index);

  } else {
    // Casting to int for MPI_Scatterv
    for (int i = 0; i < nproc; i++) {
      send_counts_int[i] = (int)send_counts[i];
      displs_int[i] = (int)displs[i];
    }

    MPI_Scatterv(graph, send_counts_int, displs_int, MPI_EDGE, *local_graph, send_counts_int[rank],
                 MPI_EDGE, 0, active_comm);

    if (*local_graph == NULL) {
      printf("Error distributing data.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  free(send_counts_int);
  free(displs_int);
}


void merge_graph(Edge** local_graph, int64_t send_counts, int nproc, int rank,
                   MPI_Comm active_comm) {
  int power = 1;
  while (nproc > 1) {
    if ((rank / power) % 2 == 0) {
      int recv_rank = rank + power;
      int recv_count = 0;

      MPI_Recv(&recv_count, 1, MPI_INT, recv_rank, 0, active_comm, MPI_STATUS_IGNORE);

      Edge* recv_graph = (Edge*)malloc(sizeof(Edge) * recv_count);
      if (recv_graph == NULL) {
        printf("Memory allocation failed for received graph.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      MPI_Recv(recv_graph, recv_count*sizeof(Edge), MPI_BYTE, recv_rank, 0, active_comm, MPI_STATUS_IGNORE);
      send_counts += recv_count;
      (*local_graph) = (Edge*)realloc(*local_graph, sizeof(Edge) * send_counts);
      if (*local_graph == NULL) {
        printf("Memory allocation failed for local graph.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // Inverted merge
      int i = send_counts - 1;
      int j = send_counts - recv_count - 1;
      int k = recv_count - 1;

      while (k >= 0) {
        if (j >= 0 && (*local_graph)[j].weight > recv_graph[k].weight) {
          (*local_graph)[i--] = (*local_graph)[j--];
        } else {
          (*local_graph)[i--] = recv_graph[k--];
        }
      }
      free(recv_graph);
    } else {
      int send_rank = rank - power;
      MPI_Send(&send_counts, 1, MPI_INT, send_rank, 0, active_comm);
      MPI_Send(*local_graph, send_counts*sizeof(Edge), MPI_BYTE, send_rank, 0, active_comm);
      break;
    }
    power *= 2;
    nproc /= 2;
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

    double distribute_unordered_percentage = (times[i][9].distribute_unordered_graph_time / total) * 100.0;
    double qsort_percentage = (times[i][9].qsort_time / total) * 100.0;
    double merge_graph_percentage = (times[i][9].merge_graph_time / total) * 100.0;
    double distribute_ordered_percentage = (times[i][9].distribute_ordered_graph_time / total) * 100.0;
    double kruskal_percentage = (times[i][9].kruskal_time / total) * 100.0;
    double collect_mst_percentage = (times[i][9].collect_mst_time / total) * 100.0;

    printf("| DUG %6.2f%% | QS %6.2f%% | MG %6.2f%% | DOG %6.2f%% | K %6.2f%% | CM %6.2f%%\n",
           distribute_unordered_percentage, qsort_percentage, merge_graph_percentage,
           distribute_ordered_percentage, kruskal_percentage, collect_mst_percentage);
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

  Edge* local_graph = NULL;
  Edge* local_mst = NULL;
  int64_t* send_counts = (int64_t*)malloc(sizeof(int64_t) * nproc);
  int64_t* displs = (int64_t*)malloc(sizeof(int64_t) * nproc);

  if (!send_counts || !displs) {
    printf("Memory allocation failed for send counts and displacements.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int result_mst_wgt = 0;
  // Calculate the number of times the algorithm will run
  // Sequential Kruskal + Parallel Kruskal with 2, 4, 8, ... processes
  int algo_times = 1 + (int)log2(nproc);
  AlgoTimes times[algo_times][10];  // 10 graph densities from 10% to 100%

  double total_start_time, total_end_time;
  double distribute_unordered_start, distribute_unordered_end;
  double sorting_start, sorting_end;
  double merge_graph_start, merge_graph_end;
  double distribute_ordered_start, distribute_ordered_end;
  double kruskal_execution_start, kruskal_execution_end;
  double collect_mst_start, collect_mst_end;

  get_argv(argc, argv, &V);
  int max_wgt = INT32_MAX / (V - 1);

  srand(time(NULL));

  // Broadcast the number of vertices
  MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Graph density from 10% to 100%
  for (int i = 1; i <= 10; i++) {    
    if (rank == 0) {
      printf("===== Density: %d%% =====\n", i * 10);
      fflush(stdout);
    }

    if (i < 10) {
      E = V - 1 + calculate_missing_edges(V, i);
      if (E <= V - 1)
        continue;  // Skip densities with less than V - 1 edges (not connected graph)

      if (rank == 0) {  // Process 0 generates the graph
        sparse_graph(&graph, V, i, max_wgt);
      }
    } else {
      E = (int64_t)V * (V - 1) / 2;
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

    // =========== Sequential Kruskal Algorithm ===========
    if (rank == 0) {
      printf("Sequential Kruskal: ...\n");
      fflush(stdout);

      sets = (Subset*)malloc(sizeof(Subset) * V);
      if (sets == NULL) {
        printf("Memory allocation failed for sets.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      total_start_time = MPI_Wtime();
      sorting_start = total_start_time;

      printf("Let's do qsort...\n");
      fflush(stdout);

      qsort(graph, E, sizeof(Edge), compare_edges);

      sorting_end = MPI_Wtime();
      kruskal_execution_start = sorting_end;

      printf("Let's do kruskal...\n");
      fflush(stdout);

      kruskal(graph, V, E, &mst, &mst_edges, &mst_weight, sets);

      kruskal_execution_end = MPI_Wtime();
      total_end_time = kruskal_execution_end;

      printf("Set result...\n");
      fflush(stdout);

      result_mst_wgt = mst_weight;
      times[0][i - 1].distribute_unordered_graph_time = 0;
      times[0][i - 1].qsort_time = sorting_end - sorting_start;
      times[0][i - 1].merge_graph_time = 0;
      times[0][i - 1].distribute_ordered_graph_time = 0;
      times[0][i - 1].kruskal_time = kruskal_execution_end - kruskal_execution_start;
      times[0][i - 1].collect_mst_time = 0;
      times[0][i - 1].total_time = total_end_time - total_start_time;

      free(sets);
      free(mst);

      printf("Done\n\n");
      fflush(stdout);
    }

    // =========== Parallel Kruskal Algorithm ===========
    // Iterate the number of processes, increasing by powers of 2 (2, 4, 8, ...)
    for (int j = 2; j <= nproc; j *= 2) {
      if (rank == 0) {
        printf("Parallel Kruskal: %d processes ...\n", j);
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
        total_start_time = MPI_Wtime();
        distribute_unordered_start = total_start_time;

        printf("Let's distribute the unordered graph...\n");
        fflush(stdout);
      }

      distibute_data(graph, &local_graph, E, send_counts, displs, j, rank, active_comm);

      if (rank == 0) {
        distribute_unordered_end = MPI_Wtime();
        sorting_start = distribute_unordered_end;

        printf("Let's do qsort...\n");
        fflush(stdout);
      }

      qsort(local_graph, send_counts[rank], sizeof(Edge), compare_edges);

      if (rank == 0) {
        sorting_end = MPI_Wtime();
        merge_graph_start = sorting_end;

        printf("Let's collect the ordered graph...\n");
        fflush(stdout);
      }

      local_graph = (Edge*)realloc(local_graph, sizeof(Edge) * send_counts[rank]);
      if (local_graph == NULL) {
        printf("Memory allocation failed for local graph.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      merge_graph(&local_graph, send_counts[rank], j, rank, active_comm);
      
      if (rank == 0) {
        memcpy(graph, local_graph, sizeof(Edge) * E);
        free(local_graph);

        merge_graph_end = MPI_Wtime();
        distribute_ordered_start = merge_graph_end;

        printf("Let's distribute the ordered graph...\n");
        fflush(stdout);
      }

      // Distribute the data to the processes
      distibute_data(graph, &local_graph, E, send_counts, displs, j, rank, active_comm);

      if (rank == 0) {
        distribute_ordered_end = MPI_Wtime();
        kruskal_execution_start = distribute_ordered_end;

        printf("Let's do kruskal...\n");
        fflush(stdout);
      }

      // Compute the local MST
      kruskal(local_graph, V, send_counts[rank], &local_mst, &mst_edges, &mst_weight, sets);

      if (rank == 0) {
        kruskal_execution_end = MPI_Wtime();
        collect_mst_start = kruskal_execution_end;

        printf("Let's collect local mst...\n");
        fflush(stdout);
      }

      // Collect the MST from all the processes
      collect_mst(rank, j, V, &local_mst, &mst_edges, &mst_weight, sets, active_comm);

      if (rank == 0) {
        collect_mst_end = MPI_Wtime();
        total_end_time = collect_mst_end;

        printf("Set result...\n");
        fflush(stdout);

        int run_index = (int)log2(j);
        times[run_index][i - 1].distribute_unordered_graph_time = distribute_unordered_end - distribute_unordered_start;
        times[run_index][i - 1].qsort_time = sorting_end - sorting_start;
        times[run_index][i - 1].merge_graph_time = merge_graph_end - merge_graph_start;
        times[run_index][i - 1].distribute_ordered_graph_time = distribute_ordered_end - distribute_ordered_start;
        times[run_index][i - 1].kruskal_time = kruskal_execution_end - kruskal_execution_start;
        times[run_index][i - 1].collect_mst_time = collect_mst_end - collect_mst_start;
        times[run_index][i - 1].total_time = total_end_time - total_start_time;

        if (mst_weight != result_mst_wgt) {
          printf("Error: MST weight is different from sequential.\n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        printf("Done\n\n");
        fflush(stdout);
      }

      free(sets);
      free(local_graph);
      free(local_mst);
    }

    if (rank == 0) {
      free(original_graph);
      free(graph);
    }
  }

  // Print the results
  if (rank == 0) {
    printf("\n\n");
    printf(" ========== Results ==========\n");

    print_field(times, algo_times, "Distribute Unordered Graph Time", get_distribute_unordered_graph_time);
    print_field(times, algo_times, "QSort Time", get_qsort_time);
    print_field(times, algo_times, "Merge Graph Time", get_merge_graph_time);
    print_field(times, algo_times, "Distribute Ordered Graph Time", get_distribute_ordered_graph_time);
    print_field(times, algo_times, "Kruskal Time", get_kruskal_time);
    print_field(times, algo_times, "Collect MST Time", get_collect_mst_time);
    print_field(times, algo_times, "Total Time", get_total_time);

    print_evaluation(times, algo_times);
  }

  MPI_Type_free(&MPI_EDGE);
  MPI_Finalize();

  return 0;
}