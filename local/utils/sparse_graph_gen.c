#include "graph_gen_utils.h"

// Add random edges to the graph
void add_random_edges(int* adj_matrix, int V, int more_edges, int max_wgt, FILE* fp) {
  int i = 0;
  while (i < more_edges) {
    int u = rand() % V;
    int v = rand() % V;
    int w = 1 + rand() % max_wgt;

    if (u != v && adj_matrix[u * V + v] == 0) {
      fprintf(fp, "%d %d %d\n", u, v, w);
      adj_matrix[u * V + v] = 1;
      adj_matrix[v * V + u] = 1;
      i++;
    }
  }
}

// Generate a connected graph
void generate_graph(int V, int more_edges, int max_wgt, FILE* fp) {
  int* nodes = (int*)malloc(sizeof(int) * V);
  int* adj_matrix = (int*)calloc(V * V, sizeof(int));

  if (!nodes || !adj_matrix) {
    printf("Memory allocation failed.\n");
    return;
  }

  // Initialize the nodes
  for (int i = 0; i < V; i++) {
    nodes[i] = i;
  }

  // Shuffle the nodes
  shuffle(nodes, V);

  fprintf(fp, "%d %d\n", V, (V - 1) + more_edges);

  // Connect all nodes
  for (int i = 1; i < V; i++) {
    int u = nodes[i];
    int v = nodes[rand() % i];
    int w = 1 + rand() % max_wgt;

    fprintf(fp, "%d %d %d\n", u, v, w);

    adj_matrix[u * V + v] = 1;
    adj_matrix[v * V + u] = 1;
  }	// All nodes are now connected

  // Add random edges
  add_random_edges(adj_matrix, V, more_edges, max_wgt, fp);

  free(adj_matrix);
  free(nodes);
}

int main(int argc, char* argv[]) {
  int V, max_wgt;
  char filename[100];
  FILE* fp;

  get_params(argc, argv, &V);

  max_wgt = INT32_MAX / (V - 1);

	srand(time(NULL));
	// Add random number of edges to the already connected graph
  int more_edges = random_num_gen(1, V * (V - 1) / 2 - (V - 1));
  double density = calculate_density(V, V - 1 + more_edges);

  // Generate filename with number of nodes and density
  sprintf(filename, "../graphs/graph_%d_nodes_%.2f_density.txt", V, density);

  open_file(&fp, filename, "w");
  generate_graph(V, more_edges, max_wgt, fp);

  printf("Graph with %d vertices and %.2f density generated.\n", V, density);

  fclose(fp);
  return 0;
}