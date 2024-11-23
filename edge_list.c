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

Edge* read_graph(int* V, int* E);
Edge* kruskal(Edge* graph, int V, int E, int* mst_weight, int* mst_edges);
int comprare_edges(const void* a, const void* b);
int find(Subset* sets, int i);
void union_sets(Subset* sets, int u, int v);

int main() {
  int V = 0, E = 0;
  int mst_weight = 0;
  int mst_edges = 0;
  Edge* graph = NULL;

  graph = read_graph(&V, &E);

  if (graph == NULL) {
    return 1;
  }

  Edge* mst = kruskal(graph, V, E, &mst_weight, &mst_edges);

  printf("MST weight: %d\n", mst_weight);
  printf("MST edges: %d\n", mst_edges);
  for (int i = 0; i < mst_edges; i++) {
    printf("%d %d %d\n", mst[i].u, mst[i].v, mst[i].weight);
  }

  free(graph);
  free(mst);
  return 0;
}

Edge* read_graph(int* V, int* E) {
  FILE* fp;
  const char* in_file = "graph.txt";
  int u, v, w;

  if ((fp = fopen(in_file, "r")) == NULL) {
    printf("Unable to open file %s for reading.\n", in_file);
    return NULL;
  }

  // Read the number of vertices
  fscanf(fp, "%d", V);
  *E = (*V) * ((*V) - 1) / 2;
  Edge* graph = (Edge*)malloc(sizeof(Edge) * (*E));

  if (graph == NULL) {
    printf("Memory allocation failed.\n");
    return NULL;
  }

  int i = 0;
  while (fscanf(fp, "%d %d %d", &u, &v, &w) != EOF) {
    graph[i].u = u;
    graph[i].v = v;
    graph[i].weight = w;
    i++;
  }

  fclose(fp);
  return graph;
}

Edge* kruskal(Edge* graph, int V, int E, int* mst_weight, int* mst_edges) {
  // Sort the edges
  qsort(graph, E, sizeof(Edge), comprare_edges);

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

int comprare_edges(const void* a, const void* b) {
  return ((Edge*)a)->weight - ((Edge*)b)->weight;
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