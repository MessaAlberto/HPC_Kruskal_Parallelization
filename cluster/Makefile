# Directories
BIN_FOLDER := bin
OUT_FOLDER := output
ERR_FOLDER := output_err
SCRIPT_FOLDER := script

# Files
MAIN_SRC = kruskal_pivotSort.c
# MAIN_SRC = kruskal_mergeSort.c
MAIN_BIN = kruskal

# Compiler and flags
CC = mpicc
CFLAGS = -Wall -Wextra -g -std=c99

# Rules
all: env $(BIN_FOLDER)/$(MAIN_BIN) test

env: create_dir
	chmod +x $(SCRIPT_FOLDER)/*.sh

create_dir:
	mkdir -p $(OUT_FOLDER) $(ERR_FOLDER) $(BIN_FOLDER)

$(BIN_FOLDER)/$(MAIN_BIN): $(MAIN_SRC)
	module load mpich-3.2 && $(CC) $(CFLAGS) -o $@ $^ -lm

test:
	$(SCRIPT_FOLDER)/test.sh

clean: cl_test
	rm -f $(BIN_FOLDER)/* || true

cl_test:
	rm -f $(OUT_FOLDER)/*.out $(ERR_FOLDER)/*.err || true