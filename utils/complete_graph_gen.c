#include "graph_gen_utils.h"

void complete_graph(int V, int max_wgt, FILE* fp) {
  fprintf(fp, "%d %d\n", V, V * (V - 1) / 2);

  for (int i = 0; i < V; i++) {
    for (int j = i + 1; j < V; j++) {
      fprintf(fp, "%d %d %d\n", i, j, 1 + rand() % max_wgt);
    }
  }
}

int main(int argc, char* argv[]) {
  int V, max_wgt;
  char filename[100];
  FILE* fp;

  get_params(argc, argv, &V, &max_wgt);

  srand(time(NULL));

  // Generate filename with number of vertices and density
  sprintf(filename, "../graphs/graph_%d_nodes_%.2f_density.txt", V, 1.0);

  open_file(&fp, filename, "w");
  complete_graph(V, max_wgt, fp);

  printf("Complete graph with %d vertices generated.\n", V);

  fclose(fp);
  return 0;
}
