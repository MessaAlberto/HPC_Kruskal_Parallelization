OUTPUT_DIR = ../output
ERROR_DIR = ../output_err

CC = mpicc
CFLAGS = -Wall -Wextra -g -std=c99

OUTPUT = kruskal
SRC = kruskal.c
OBJS = $(SRC:.c=.o)

all: env test

env: create_output_dir module_load $(OUTPUT)
	chmod +x test.sh qsub_dir.sh

create_output_dir:
	mkdir -p $(OUTPUT_DIR) $(ERROR_DIR)

module_load:
	$(SHELL) -l -c "module load mpich-3.2"

$(OUTPUT): $(OBJS)
	$(CC) $(CFLAGS) -o $(OUTPUT) $(OBJS) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: 
	./test.sh

clean: cl_test
	rm -f $(OUTPUT) $(OBJS)

cl_test:
	rm -f $(OUTPUT_DIR)/*.out $(ERROR_DIR)/*.err || true