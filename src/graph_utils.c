#include "graph_utils.h"

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

int compare_edges(const void* a, const void* b) {
  return ((Edge*)a)->weight - ((Edge*)b)->weight;
}

void create_output_filename(char* input_filename, char* prefix, char** output_filename) { 
  char* base_name = basename(input_filename); 

  *output_filename = (char*)malloc(strlen(base_name) + strlen(prefix) + strlen("../mst/") + 5);
  if (*output_filename == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  sprintf(*output_filename, "../mst/%s%s", prefix, base_name);
}