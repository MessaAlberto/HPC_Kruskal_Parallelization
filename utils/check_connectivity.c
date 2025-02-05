#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// DFS
void DFS(int* adj_matrix, int v, int* visited, int node) {
  visited[node] = 1;  // Mark the node as visited

  // Visit all the adjacent nodes
  for (int i = 0; i < v; i++) {
    if (*(adj_matrix + node * v + i) == 1 && !visited[i]) {
      DFS(adj_matrix, v, visited, i);
    }
  }
}

// Check if the graph is connected
int is_connected(int* adj_matrix, int v) {
  int* visited = (int*)malloc(sizeof(int) * v);

  if (visited == NULL) {
    printf("Memory allocation failed!\n");
    return 0;
  }

  // Initialize the visited array with 0
  for (int i = 0; i < v; i++) {
    visited[i] = 0;
  }

  DFS(adj_matrix, v, visited, 0);

  // Check if all the nodes are visited
  for (int i = 0; i < v; i++) {
    if (visited[i] == 0) {
      free(visited);
      return 0;
    }
  }

  // If all the nodes are visited, the graph is connected
  free(visited);
  return 1;
}

int main(int argc, char* argv[]) {
  int V, E;
  FILE* fp;

  if (argc < 2) {
    printf("Usage: %s ../graphs/<input_file>\n", argv[0]);
    return 1;
  }

  char* filename = argv[1];

  if ((fp = fopen(filename, "r")) == NULL) {
    printf("Unable to open file %s for reading.\n", filename);
    return 1;
  }

  // Read the number of vertices and edges
  fscanf(fp, "%d %d", &V, &E);

  // Allocate memory for the adjacency matrix
  int* adj_matrix = (int*)malloc(sizeof(int) * V * V);

  if (adj_matrix == NULL) {
    printf("Memory allocation failed!\n");
    return 1;
  }

  // Initialize the adjacency matrix with 0
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < V; j++) {
      *(adj_matrix + i * V + j) = 0;
    }
  }

  // Read the edges from the file
  int u, w, z;
  for (int i = 0; i < E; i++) {
    fscanf(fp, "%d %d %d", &u, &w, &z);
    *(adj_matrix + u * V + w) = 1;
    *(adj_matrix + w * V + u) = 1;
  }

  if (is_connected(adj_matrix, V)) {
    printf("The graph is connected.\n");
  } else {
    printf("The graph is not connected.\n");
  }

  free(adj_matrix);

  return 0;
}
