#include <mpi.h>

#include "disjoint_set.h"
#include "graph_utils.h"

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

void read_input(int argc, char* argv[], Edge** graph, int* V, int* E,
                FILE** fp) {
  if (argc < 2) {
    printf("Usage: %s <input_file>\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  *graph = read_graph(V, E, argv[1]);

  if (*graph == NULL) {
    printf("Error reading graph.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  *fp = fopen("output.txt", "w");
  if (*fp == NULL) {
    printf("Error opening output file.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void distibute_data(int E, int nproc, int* sendCounts, int* displs, Edge* graph,
                    Edge** local_graph, int rank) {
  int base_size = E / nproc;
  int extra = E % nproc;

  int offset = 0;
  for (int i = 0; i < nproc; i++) {
    sendCounts[i] = base_size + (i < extra ? 1 : 0);
    displs[i] = offset;
    offset += sendCounts[i];
  }

  int local_E = sendCounts[rank];
  *local_graph = (Edge*)malloc(sizeof(Edge) * local_E);
  if (*local_graph == NULL) {
    printf("Memory allocation failed for local graph.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Scatterv(graph, sendCounts, displs, MPI_EDGE, *local_graph, local_E,
               MPI_EDGE, 0, MPI_COMM_WORLD);

  if (local_graph == NULL) {
    printf("Error distributing data.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

Edge* kruskal(Edge* graph, int V, int E, int* mst_weight, int* mst_edges,
              Subset* sets) {
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

  for (int i = 0; i < E && *mst_edges < V - 1; i++) {
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

void kruskal_mst(Edge* local_graph, int V, int local_E, int* mst_weight,
                 int* mst_edges, Subset* sets, Edge** local_mst) {
  *local_mst = kruskal(local_graph, V, local_E, mst_weight, mst_edges, sets);
  if (*local_mst == NULL) {
    printf("Error computing MST.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void collect_mst(int rank, int nproc, int mst_edges, Edge** local_mst,
                 int mst_weight, Subset* sets) {
  int power = 1;
  while (nproc > 1) {
    if ((rank / power) % 2 == 0) {
      int recv_rank = rank + power;
      if (recv_rank < nproc) {
        int recv_E = 0;
        MPI_Recv(&recv_E, 1, MPI_INT, recv_rank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        Edge* recv_mst = (Edge*)malloc(sizeof(Edge) * recv_E);
        if (recv_mst == NULL) {
          printf("Memory allocation failed for received MST.\n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        MPI_Recv(recv_mst, recv_E, MPI_EDGE, recv_rank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        Edge* new_mst = (Edge*)malloc(sizeof(Edge) * (mst_edges + recv_E));
        if (new_mst == NULL) {
          printf("Memory allocation failed for new MST.\n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int i = 0;
        for (i = 0; i < mst_edges; i++) {
          new_mst[i] = (*local_mst)[i];
        }

        int count_new_edges = 0;
        for (int j = 0; j < recv_E; j++) {
          int u = recv_mst[j].u;
          int v = recv_mst[j].v;

          int root_u = find(sets, u);
          int root_v = find(sets, v);

          if (root_u != root_v) {
            new_mst[i++] = recv_mst[j];
            mst_weight += recv_mst[j].weight;
            union_sets(sets, root_u, root_v);
            count_new_edges++;
          }
        }

        mst_edges += count_new_edges;
        free(*local_mst);
        free(recv_mst);
        *local_mst = new_mst;
      }
    } else {
      int send_rank = rank - power;
      MPI_Send(&mst_edges, 1, MPI_INT, send_rank, 0, MPI_COMM_WORLD);
      MPI_Send(*local_mst, mst_edges, MPI_EDGE, send_rank, 0, MPI_COMM_WORLD);
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

  int V = 0, E = 0;
  int mst_weight = 0, mst_edges = 0;
  Edge* graph = NULL;
  FILE* fp;
  char* output_filename;

  if (rank == 0) {
    read_input(argc, argv, &graph, &V, &E, &fp);
    create_output_filename(argv[1], "mpi_mst_", &output_filename);
    printf("output_filename: %s\n", output_filename);
    qsort(graph, E, sizeof(Edge), compare_edges);
  }

  MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&E, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int sendCounts[nproc], displs[nproc];
  Edge* local_graph;

  distibute_data(E, nproc, sendCounts, displs, graph, &local_graph, rank);

  int local_E = sendCounts[rank];
  printf("Rank: %d, local_E: %d\n", rank, local_E);

  if (rank == 0) free(graph);

  Subset* sets = (Subset*)malloc(sizeof(Subset) * V);
  if (sets == NULL) {
    printf("Memory allocation failed for sets.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  Edge* local_mst = NULL;
  kruskal_mst(local_graph, V, sendCounts[rank], &mst_weight, &mst_edges, sets,
              &local_mst);

  collect_mst(rank, nproc, mst_edges, &local_mst, mst_weight, sets);

  // print mst
  if (rank == 0) {
    FILE* out_fp = fopen(output_filename, "w");
    if (out_fp == NULL) {
      printf("Error opening output file.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fprintf(out_fp, "%d %d %d\n", V, mst_edges, mst_weight);
    for (int i = 0; i < mst_edges; i++) {
      fprintf(out_fp, "%d %d %d\n", local_mst[i].u, local_mst[i].v,
              local_mst[i].weight);
    }
  }

  free(local_graph);
  free(sets);
  free(local_mst);

  if (rank == 0) fclose(fp);

  MPI_Type_free(&MPI_EDGE);
  MPI_Finalize();

  return 0;
}