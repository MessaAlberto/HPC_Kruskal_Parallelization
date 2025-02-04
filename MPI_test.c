#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Edge {
  int u;
  int v;
  int weight;
} Edge;

typedef struct Subset {
  int parent;
  int rank;
} Subset;

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

int comprare_edges(const void* a, const void* b) {
  return ((Edge*)a)->weight - ((Edge*)b)->weight;
}

Edge* read_graph(int* V, int* E, char* in_file) {
  FILE* fp;
  int u, v, w;

  if ((fp = fopen(in_file, "r")) == NULL) {
    printf("Unable to open file %s for reading.\n", in_file);
    return NULL;
  }

  if (fscanf(fp, "%d %d", V, E) != 2) {
    printf("Invalid input file format.\n");
    fclose(fp);
    return NULL;
  }

  Edge* graph = (Edge*)malloc(sizeof(Edge) * (*E));

  if (graph == NULL) {
    printf("Memory allocation failed.\n");
    fclose(fp);
    return NULL;
  }

  for (int i = 0; i < *E; i++) {
    if (fscanf(fp, "%d %d %d", &u, &v, &w) != 3) {
      printf("Invalid input file format.\n");
      free(graph);
      fclose(fp);
      return NULL;
    }

    graph[i].u = u;
    graph[i].v = v;
    graph[i].weight = w;
  }

  fclose(fp);
  return graph;
}

int find(Subset* sets, int i) {
  if (sets[i].parent != i) {
    sets[i].parent = find(sets, sets[i].parent);
  }
  return sets[i].parent;
}

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

Edge* kruskal(Edge* graph, int V, int E, int* mst_weight, int* mst_edges, Subset* sets) {
  if (sets == NULL) {
    printf("Memory allocation failed.\n");
    return NULL;
  }

  for (int i = 0; i < V; i++) {
    sets[i].parent = i;
    sets[i].rank = 0;
  }

  // // MST data structure
  Edge* mst = (Edge*)malloc(sizeof(Edge) * (V - 1));
  *mst_weight = 0;
  *mst_edges = 0;

  for (int i = 0; i < E && *mst_edges < V - 1; i++) {
    int u = graph[i].u;
    int v = graph[i].v;

    int root_u = find(sets, u);
    int root_v = find(sets, v);

    if (root_u != root_v) {
      // Include the edge in the MST
      mst[(*mst_edges)++] = graph[i];
      *mst_weight += graph[i].weight;

      // Union the sets
      union_sets(sets, root_u, root_v);
    }
  }
  return mst;
}

int main(int argc, char* argv[]) {
  // MPI init
  MPI_Init(&argc, &argv);

  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  create_mpi_edge_type();

  // Get input
  Edge* graph = NULL;
  int V = 0, E = 0;
  int mst_weight = 0;
  int mst_edges = 0;

  // create ouput file
  FILE* fp;
  if ((fp = fopen("output.txt", "w")) == NULL) {
    printf("Unable to open file %s for writing.\n", "output.txt");
    return 1;
  }
  // print my rank
  printf("Rank: %d\n", rank);

  if (rank == 0) {
    if (argc < 2) {
      printf("Usage: %s <input_file>\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    graph = read_graph(&V, &E, argv[1]);

    if (graph == NULL) {
      printf("Error reading graph.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fprintf(fp, "%d %d\n", V, E);

    // // Quick sort
    qsort(graph, E, sizeof(Edge), comprare_edges);

    // fprintf(fp, "Sorted graph:\n");
    // for (int i = 0; i < E; i++) {
    //   fprintf(fp, "%d %d %d\n", graph[i].u, graph[i].v, graph[i].weight);
    // }
  }

  MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&E, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // // Parsing sorted graph to processes
  int sendCounts[nproc];
  int displs[nproc];

  int base_size = E / nproc;
  int extra = E % nproc;

  int offset = 0;
  for (int i = 0; i < nproc; i++) {
    sendCounts[i] = base_size + (i < extra ? 1 : 0);
    displs[i] = offset;
    offset += sendCounts[i];
  }
  
  int local_E = sendCounts[rank];
  Edge* local_graph = (Edge*)malloc(sizeof(Edge) * local_E);

  MPI_Scatterv(graph, sendCounts, displs, MPI_EDGE, local_graph, local_E, MPI_EDGE, 0, MPI_COMM_WORLD);

  printf("Rank: %d, local_E: %d\n", rank, local_E);

  if (rank == 0) free(graph);
  

  // kruskal
  // Union-Find data structure
  Subset* sets = (Subset*)malloc(sizeof(Subset) * V);
  Edge* local_mst = (Edge*)malloc(sizeof(Edge) * local_E);
  local_mst = kruskal(local_graph, V, local_E, &mst_weight, &mst_edges, sets);


  // while until 1 process {
  // If proc mod something == 0
  //   send
  // else
  //   receive
  //   mstMerge
  // }

  int power = 1;
  while (nproc > 1) {
    
    if ((rank / power) % 2 == 0) {
      int recv_rank = rank + power;
      if (recv_rank < nproc) {
        int recv_E = 0;
        MPI_Recv(&recv_E, 1, MPI_INT, recv_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        Edge* recv_mst = (Edge*)malloc(sizeof(Edge) * recv_E);
        MPI_Recv(recv_mst, recv_E, MPI_EDGE, recv_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Merge mst
        Edge* new_mst = (Edge*)malloc(sizeof(Edge) * (mst_edges + recv_E));
        int i = 0;
        for (i = 0; i < mst_edges; i++) {
          new_mst[i] = local_mst[i];
        }

        int count_new_edges = 0;
        for (int j = 0; j < recv_E; j++) {
          // if edge connect two different sets then add it to mst
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
        free(local_mst);
        local_mst = new_mst;
      }
    } else {
      int send_rank = rank - power;
      MPI_Send(&mst_edges, 1, MPI_INT, send_rank, 0, MPI_COMM_WORLD);
      MPI_Send(local_mst, mst_edges, MPI_EDGE, send_rank, 0, MPI_COMM_WORLD);
      free(local_mst);
      free(sets);
      break;
    }
    power *= 2;
    nproc /= 2;
  } // one process remaining

  // print mst
  if (rank == 0) {
    fprintf(fp, "MST:\n");
    for (int i = 0; i < mst_edges; i++) {
      fprintf(fp, "%d %d %d\n", local_mst[i].u, local_mst[i].v, local_mst[i].weight);
    }
    fprintf(fp, "MST weight: %d\n", mst_weight);
  }
  
  free(local_graph);

  MPI_Type_free(&MPI_EDGE);
  MPI_Finalize();
}