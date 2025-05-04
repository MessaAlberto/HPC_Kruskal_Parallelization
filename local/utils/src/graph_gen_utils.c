#include "graph_gen_utils.h"

void open_file(FILE** fp, const char* filename, const char* mode) {
  if ((*fp = fopen(filename, mode)) == NULL) {
    printf("Unable to open file %s for writing.\n", filename);
    exit(1);
  }
}

// Get the number of vertices and maximum weight from the command line
void get_params(int argc, char* argv[], int* V) {
  if (argc != 2) {
    printf("Usage: ./graph_generator <number_of_vertices>\n");
    exit(1);
  }

  *V = atoi(argv[1]);

  if (*V <= 0) {
    printf("Number of vertices must be greater than 0.\n");
    exit(1);
  }
}

// Generate a random number between min and max (inclusive)
long long random_num_gen(long long min, long long max) {
    if (min > max) {
        printf("Error: min (%lld) is greater than max (%lld)\n", min, max);
        exit(1);
    }

    long long range = max - min + 1;
    long long rand_val;

    do {
        rand_val = ((long long)rand() << 16) + rand();
    } while (rand_val >= (RAND_MAX * RAND_MAX) / range * range);

    return min + (rand_val % range);
}

// Shuffle the nodes
void shuffle(int* nodes, int V) {
  for (int i = 0; i < V; i++) {
    int j = rand() % V;
    int temp = nodes[i];
    nodes[i] = nodes[j];
    nodes[j] = temp;
  }
}

double calculate_density(int V, int E) {
  return (2.0 * E) / (V * (V - 1));
}