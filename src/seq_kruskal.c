#include "disjoint_set.h"
#include "graph_utils.h"

void read_input(int argc, char* argv[], Edge** graph, int* V, int* E, FILE** fp) {
  if (argc < 2) {
    printf("Usage: %s <input_file>\n", argv[0]);
    exit(1);
  }

  *graph = read_graph(V, E, argv[1]);

  if (*graph == NULL) {
    printf("Error reading graph.\n");
    exit(1);
  }

  *fp = fopen("output.txt", "w");
  if (*fp == NULL) {
    printf("Error opening output file.\n");
    exit(1);
  }
}

Edge* kruskal(Edge* graph, int V, int E, int* mst_weight, int* mst_edges) {
  // Sort the edges
  qsort(graph, E, sizeof(Edge), compare_edges);

  // Union-Find data structure
  Subset* sets = (Subset*)malloc(sizeof(Subset) * V);

  if (sets == NULL) {
    printf("Memory allocation failed.\n");
    return NULL;
  }

  for (int i = 0; i < V; i++) {
    sets[i].parent = i;
    sets[i].rank = 0;
  }

  // MST data structure
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

  free(sets); 
  return mst;
}

int main(int argc, char* argv[]) {
  int V = 0, E = 0;
  int mst_weight = 0, mst_edges = 0;
  Edge* graph = NULL;
  FILE* fp;
  char* output_filename;

  read_input(argc, argv, &graph, &V, &E, &fp);
  create_output_filename(argv[1], "seq_mst_", &output_filename);

  Edge* mst = kruskal(graph, V, E, &mst_weight, &mst_edges);

  // Output file
  FILE* out_fp = fopen(output_filename, "w");
  if (out_fp == NULL) {
    printf("Error opening output file.\n");
    exit(1);
  }

  fprintf(out_fp, "%d %d %d\n", V, mst_edges, mst_weight);
  for (int i = 0; i < mst_edges; i++) {
    fprintf(out_fp, "%d %d %d\n", mst[i].u, mst[i].v, mst[i].weight);
  }

  fclose(out_fp);
  fclose(fp);
  free(graph);
  free(mst);
  free(output_filename);
  return 0;
}