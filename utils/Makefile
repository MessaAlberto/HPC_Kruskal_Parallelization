# Directory for generated graphs
GRAPH_DIR = ../graphs

# Compiler and flags
CC = gcc
CFLAGS = -Wall -O2

# Source files
SRC_GRAPH_GEN = graph_gen.c graph_gen_utils.c
SRC_COMPLETE_GRAPH = complete_graph_gen.c graph_gen_utils.c
SRC_CHECK_CONNECTIVITY = check_connectivity.c

# Object files
OBJ_GRAPH_GEN = $(SRC_GRAPH_GEN:.c=.o)
OBJ_COMPLETE_GRAPH = $(SRC_COMPLETE_GRAPH:.c=.o)
OBJ_CHECK_CONNECTIVITY = $(SRC_CHECK_CONNECTIVITY:.c=.o)

# Executables
BIN_GRAPH_GEN = graph_gen.exe
BIN_COMPLETE_GRAPH = complete_graph_gen.exe
BIN_CHECK_CONNECTIVITY = check_connectivity.exe

# Default target: build all
all: $(GRAPH_DIR) $(BIN_GRAPH_GEN) $(BIN_COMPLETE_GRAPH) $(BIN_CHECK_CONNECTIVITY)

# Rule for creating directory for generated graphs
$(GRAPH_DIR):
	mkdir -p $(GRAPH_DIR)

# Compile graph_gen
$(BIN_GRAPH_GEN): $(OBJ_GRAPH_GEN)
	$(CC) $(CFLAGS) -o $@ $^

# Compile complete_graph_gen
$(BIN_COMPLETE_GRAPH): $(OBJ_COMPLETE_GRAPH)
	$(CC) $(CFLAGS) -o $@ $^

# Compile check_connectivity
$(BIN_CHECK_CONNECTIVITY): $(OBJ_CHECK_CONNECTIVITY)
	$(CC) $(CFLAGS) -o $@ $^

# Compile C files into object files
%.o: %.c graph_gen_utils.h
	$(CC) $(CFLAGS) -c $< -o $@

# Clean compiled files
clean:
	rm -f $(BIN_GRAPH_GEN) $(BIN_COMPLETE_GRAPH) $(BIN_CHECK_CONNECTIVITY) *.o

# Clean all generated graphs
clean_graphs:
	rm -f ../graphs/*.txt

# Clean all
clean_all: clean clean_graphs