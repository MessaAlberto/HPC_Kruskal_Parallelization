# Directories
BIN_FOLDER := bin
OBJ_FOLDER := obj
SRC_FOLDER := src
INC_FOLDER := include
MST_FOLDER := ../mst

# Compiler, flags, and linker
CC := gcc
CFLAGS := -Wall -Wextra -g -I$(INC_FOLDER)
MPIFLAGS := -I "C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
MPILIB := -L "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi

# Files
SEQ_MAIN := seq_kruskal.c
PAR_MAIN := par_kruskal.c
SEQ_BIN := $(BIN_FOLDER)/seq_kruskal
PAR_BIN := $(BIN_FOLDER)/par_kruskal
COMMON_SRCS := $(wildcard $(SRC_FOLDER)/*.c)
COMMON_OBJS := $(patsubst $(SRC_FOLDER)/%.c,$(OBJ_FOLDER)/%.o,$(COMMON_SRCS))

# Rules
all: create_dir $(SEQ_BIN) $(PAR_BIN)

create_dir:
	mkdir -p $(OBJ_FOLDER) $(BIN_FOLDER) $(MST_FOLDER)

$(SEQ_BIN): $(COMMON_OBJS) $(SEQ_MAIN)
	$(CC) $(CFLAGS) $^ -o $@

$(PAR_BIN): $(COMMON_OBJS) $(PAR_MAIN)
	$(CC) $(CFLAGS) $(MPIFLAGS) $^ -o $@ $(MPILIB)

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test:
ifndef GRAPH
	$(error GRAPH is not defined. Usage: make test GRAPH=../graphs/<graph_name.txt>)
endif
	@echo "Running sequential Kruskal on $(GRAPH)"
	$(SEQ_BIN).exe $(GRAPH)
	@echo ""
	@echo "Running parallel Kruskal on $(GRAPH)"
	mpiexec -n 4 $(PAR_BIN).exe $(GRAPH)

clean:
	rm -f $(OBJ_FOLDER)/*.o $(BIN_FOLDER)/* || true

clean_mst:
	rm -f $(MST_FOLDER)/*.txt || true

clean_all: clean clean_mst