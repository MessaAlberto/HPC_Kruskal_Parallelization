# Directories
BIN_FOLDER := bin
OBJ_FOLDER := obj
SRC_FOLDER := src
INC_FOLDER := include
GRAPH_FOLDER := ../graphs

# Compiler and flags
CC := gcc
CFLAGS := -Wall -O2 -I$(INC_FOLDER)

# Files
MAINS := sparse_graph_gen complete_graph_gen check_connectivity
BINS := $(patsubst %, $(BIN_FOLDER)/%, $(MAINS))
UTIL_SRCS := $(SRC_FOLDER)/graph_gen_utils.c
UTIL_OBJS := $(OBJ_FOLDER)/graph_gen_utils.o

# Rules
all: create_dir $(BINS)

create_dir:
	mkdir -p $(BIN_FOLDER) $(OBJ_FOLDER) $(GRAPH_FOLDER)

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_FOLDER)/%: %.c $(UTIL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(BIN_FOLDER)/* $(OBJ_FOLDER)/* || true

clean_graph:
	rm -f $(GRAPH_FOLDER)/*.txt || true

clean_all: clean clean_graph