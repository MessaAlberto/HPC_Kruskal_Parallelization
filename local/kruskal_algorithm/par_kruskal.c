#include <mpi.h>

#include "disjoint_set.h"
#include "graph_utils.h"

#define MAX_WGT INT32_MAX
MPI_Datatype MPI_EDGE;


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

void read_input(int argc, char* argv[], Edge** graph, int* V, int* E) {
  if (argc < 2) {
    printf("Usage: %s <input_file>\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  *graph = read_graph(V, E, argv[1]);

  if (*graph == NULL) {
    printf("Error reading graph.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

Edge* kruskal(Edge* graph, int V, int E, int* mst_weight, int* mst_edges, Subset* sets) {
  if (sets == NULL) {
    printf("Memory allocation failed.\n");
    return NULL;
  }

  for (int i = 0; i < V; i++) {
    sets[i].parent = i;
    sets[i].rank = 0;
  }

  Edge* mst = (Edge*)malloc(sizeof(Edge) * (V - 1));
  if (mst == NULL) {
    printf("Memory allocation failed for MST.\n");
    return NULL;
  }
  *mst_weight = 0;
  *mst_edges = 0;

  for (int i = 0; *mst_edges < V - 1 && i < E; i++) {
    int u = graph[i].u;
    int v = graph[i].v;

    int root_u = find(sets, u);
    int root_v = find(sets, v);

    if (root_u != root_v) {
      mst[(*mst_edges)++] = graph[i];
      *mst_weight += graph[i].weight;
      union_sets(sets, root_u, root_v);
    }
  }

  return mst;
}

void kruskal_mst(Edge* local_graph, int V, int local_E, int* mst_weight, int* mst_edges,
                 Subset* sets, Edge** local_mst) {
  *local_mst = kruskal(local_graph, V, local_E, mst_weight, mst_edges, sets);
  if (*local_mst == NULL) {
    printf("Error computing MST.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

int pivot_sort(Edge* graph, Edge** local_graph, int E, int nproc, int rank,
               MPI_Comm active_comm, int max_wgt) {
  int local_edge_count = 0;
  int* counts = NULL;

  if (rank == 0) {
    Edge** buckets = malloc(nproc * sizeof(Edge*));
    counts = calloc(nproc, sizeof(int));
    int* bucket_caps = calloc(nproc, sizeof(int));
    int bucket_size = (E / nproc) + 1;

    for (int i = 0; i < nproc; i++) {
      buckets[i] = malloc(bucket_size * sizeof(Edge));
      bucket_caps[i] = bucket_size;
      counts[i] = 0;
    }

    int bucket_width = max_wgt / nproc;
    if (bucket_width == 0) bucket_width = 1;

    for (int i = 0; i < E; i++) {
      int bucket = graph[i].weight / bucket_width;
      if (bucket >= nproc) bucket = nproc - 1;

      if (counts[bucket] >= bucket_caps[bucket]) {
        bucket_caps[bucket] *= 2;
        buckets[bucket] = realloc(buckets[bucket], bucket_caps[bucket] * sizeof(Edge));
      }
      buckets[bucket][counts[bucket]++] = graph[i];
    }

    for (int i = 1; i < nproc; i++) {
      int count = counts[i];
      MPI_Send(&count, 1, MPI_LONG_LONG_INT, i, 0, active_comm);
      MPI_Send(buckets[i], count, MPI_EDGE, i, 1, active_comm);
    }

    local_edge_count = counts[0];
    *local_graph = malloc(local_edge_count * sizeof(Edge));
    memcpy(*local_graph, buckets[0], local_edge_count * sizeof(Edge));

    for (int i = 0; i < nproc; i++) free(buckets[i]);
    free(buckets);
    free(bucket_caps);
  } else {
    MPI_Recv(&local_edge_count, 1, MPI_LONG_LONG_INT, 0, 0, active_comm, MPI_STATUS_IGNORE);
    *local_graph = malloc(local_edge_count * sizeof(Edge));
    if (*local_graph == NULL) {
      printf("Memory allocation failed for local graph.\n");
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Recv(*local_graph, local_edge_count, MPI_EDGE, 0, 1, active_comm, MPI_STATUS_IGNORE);

    if (*local_graph == NULL) {
      printf("Error receiving local graph.\n");
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  qsort(*local_graph, local_edge_count, sizeof(Edge), compare_edges);
  if (rank == 0) {
    free(counts);
  }
  return local_edge_count;
}

void collect_mst(int rank, int nproc, int V, int* mst_edges, Edge** local_mst, int* mst_weight,
                 Subset* sets) {
  int power = 1;
  while (nproc > 1) {
    if ((rank / power) % 2 == 0) {
      int recv_rank = rank + power;

      Edge* recv_mst = (Edge*)malloc(sizeof(Edge) * (*mst_edges));

      if (recv_mst == NULL) {
        printf("Memory allocation failed for received MST.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      MPI_Recv(recv_mst, *mst_edges, MPI_EDGE, recv_rank, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      int recv_edges = *mst_edges;
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
      MPI_Send(*local_mst, *mst_edges, MPI_EDGE, send_rank, 0, MPI_COMM_WORLD);
      break;
    }
    power *= 2;
    nproc /= 2;
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  create_mpi_edge_type();

  double pivotSort_start_time, pivotSort_end_time;
  double alg_start_time, alg_end_time;
  double krustal_start_time, krustal_end_time;
  double collect_start_time, collect_end_time;

  int V = 0, E = 0;
  int mst_weight = 0, mst_edges = 0;
  Edge* graph = NULL;
  Edge* local_graph;
  int sendCounts[nproc];
  Subset* sets;
  Edge* local_mst = NULL;

  char* output_filename;

  if (rank == 0) {
    read_input(argc, argv, &graph, &V, &E);
    create_output_filename(argv[1], "mpi_mst_", &output_filename);
  }

  // Broadcast the number of vertices and edges
  MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&E, 1, MPI_INT, 0, MPI_COMM_WORLD);

  sets = (Subset*)malloc(sizeof(Subset) * V);
  if (sets == NULL) {
    printf("Memory allocation failed for sets.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Start the algorithm doing the qsort
  if (rank == 0) {
    alg_start_time = MPI_Wtime();
    pivotSort_start_time = alg_start_time;
  }

  sendCounts[rank] = pivot_sort(graph, &local_graph, E, nproc, rank, MPI_COMM_WORLD, MAX_WGT);

  if (rank == 0) {
    pivotSort_end_time = MPI_Wtime();
    krustal_start_time = pivotSort_end_time;
    free(graph);
  }

  // Compute the local MST
  kruskal_mst(local_graph, V, sendCounts[rank], &mst_weight, &mst_edges, sets, &local_mst);

  if (rank == 0) {
    krustal_end_time = MPI_Wtime();
    collect_start_time = krustal_end_time;
  }

  // Collect the MST from all the processes
  collect_mst(rank, nproc, V, &mst_edges, &local_mst, &mst_weight, sets);

  if (rank == 0) {
    collect_end_time = MPI_Wtime();
    alg_end_time = collect_end_time;

    // Write the output
    FILE* out_fp = fopen(output_filename, "w");
    if (out_fp == NULL) {
      printf("Error opening output file.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fprintf(out_fp, "%d %d %d\n", V, mst_edges, mst_weight);
    for (int i = 0; i < mst_edges; i++) {
      fprintf(out_fp, "%d %d %d\n", local_mst[i].u, local_mst[i].v, local_mst[i].weight);
    }

    fclose(out_fp);

    // Print time information
    printf("PivotSort time: %.6f seconds\n", pivotSort_end_time - pivotSort_start_time);
    printf("Kruskal time: %.6f seconds\n", krustal_end_time - krustal_start_time);
    printf("Collect time: %.6f seconds\n", collect_end_time - collect_start_time);
    printf("Tot algorithm time: %.6f seconds\n", alg_end_time - alg_start_time);
  }

  free(local_graph);
  free(sets);
  free(local_mst);

  MPI_Type_free(&MPI_EDGE);
  MPI_Finalize();

  return 0;
}