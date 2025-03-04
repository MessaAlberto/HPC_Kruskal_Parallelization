OUTPUT_DIR = ../output
ERROR_DIR = ../output_err

CC = mpicc
CFLAGS = -Wall -Wextra -g -std=c99

OUTPUT = kruskal
SRC = kruskal.c
OBJS = $(SRC:.c=.o)

all: env test

env: create_output_dir $(OUTPUT)
	chmod +x test.sh qsub_dir.sh

create_output_dir:
	mkdir -p $(OUTPUT_DIR) $(ERROR_DIR)

$(OUTPUT): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: 
	./test.sh

clean: cl_test
	rm -f $(OUTPUT) $(OBJS)

cl_test:
	rm -f $(OUTPUT_DIR)/*.out $(ERROR_DIR)/*.err || true