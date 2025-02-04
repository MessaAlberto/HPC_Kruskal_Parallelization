#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void open_file(FILE** fp, const char* filename, const char* mode) {
  if ((*fp = fopen(filename, mode)) == NULL) {
    printf("Unable to open file %s for writing.\n", filename);
    exit(1);
  }
}

void get_params(int* V, int* max_wgt) {
  printf("Number of vertices: ");
  scanf("%d", V);

  if (*V <= 0) {
    printf("Number of vertices must be greater than 0.\n");
    exit(1);
  }

  printf("Maximum weight: ");
  scanf("%d", max_wgt);

  if (*max_wgt <= 0) {
    printf("Maximum weight must be greater than 0.\n");
    exit(1);
  }
}

long long random_num_gen(long long min, long long max) {
    long long range = (long long)(max - min + 1);
    long long rand_val;

    do {
        rand_val = ((long long)rand() << 16) + rand();
    } while (rand_val >= (RAND_MAX * RAND_MAX) / range * range);

    return (long long)(min + (rand_val % range));
}

void shuffle(int* nodes, int v) {
  for (int i = 0; i < v; i++) {
    int j = rand() % v;
    int temp = nodes[i];
    nodes[i] = nodes[j];
    nodes[j] = temp;
  }
}

double calculate_connecteness(int V, int more_edges) {
  return ((double)(V - 1) + more_edges) / (V * (V - 1) / 2) * 100;
}

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

void generate_graph(int V, int max_wgt, FILE* fp) {
  int* nodes = (int*)malloc(sizeof(int) * V);
  int* adj_matrix = (int*)calloc(V * V, sizeof(int));

	// Add random number of edges
	srand(time(NULL));
	int more_edges = random_num_gen(V - 1, V * (V - 1) / 2);
	printf("Total edges: %d\n", V - 1 + more_edges);
	
	// Percentage of connecteness
  printf("Connecteness: %.2f%%\n", calculate_connecteness(V, more_edges));

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
  fprintf(fp, "%d %d %d\n", nodes[0], nodes[1], 1 + rand() % max_wgt);

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

int main() {
  int V, max_wgt;
  const char* filename = "../graph.txt";
  FILE* fp;

  open_file(&fp, filename, "w");
  get_params(&V, &max_wgt);

  generate_graph(V, max_wgt, fp);

  fclose(fp);
  return 0;
}