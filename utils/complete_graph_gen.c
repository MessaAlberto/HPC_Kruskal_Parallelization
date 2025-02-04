#include <stdio.h>
#include <stdlib.h>

void complete_graph(int v, int max_wgt, const char* out_file);

int main() {
  int v, max_wgt;
  const char* filename = "../graph.txt";

  printf("Number of vertices: ");
  scanf("%d", &v);

  if (v <= 0) {
    printf("Number of vertices must be greater than 0.\n");
    return 1;
  }

  printf("Maximum weight: ");
  scanf("%d", &max_wgt);

  if (max_wgt <= 0) {
    printf("Maximum weight must be greater than 0.\n");
    return 1;
  }

  complete_graph(v, max_wgt, filename);

  return 0;
}

void complete_graph(int v, int max_wgt, const char* out_file) {
  FILE* fp;

  if ((fp = fopen(out_file, "w")) == NULL) {
    printf("Unable to open file %s for writing.\n", out_file);
    return;
  }

  fprintf(fp, "%d\n", v);

  for (int i = 0; i < v; i++) {
    for (int j = i + 1; j < v; j++) {
      fprintf(fp, "%d %d %d\n", i, j, 1 + rand() % max_wgt);
    }
  }

  fclose(fp);

  printf("\tGraph is written to file %s.\n", out_file);
}
