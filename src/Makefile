# Directory for generated mst
MST_DIR = ../mst

# Variables for compilation and linking
CC = gcc
CFLAGS = -Wall -Wextra -g
INCLUDE = -I "C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
LIBDIR = -L "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64"
LIBS += -lmsmpi

# Source files
SEQ_SRCS = seq_kruskal.c disjoint_set.c graph_utils.c
PAR_SRCS = par_kruskal.c disjoint_set.c graph_utils.c

# Object files
SEQ_OBJS = $(SEQ_SRCS:.c=.o)
PAR_OBJS = $(PAR_SRCS:.c=.o)

# Executable file
SEQ_EXEC = seq_kruskal.exe
PAR_EXEC = par_kruskal.exe

# Rule for generating all executables
all: $(MST_DIR) $(SEQ_EXEC) $(PAR_EXEC)

# Rule for creating directory for generated mst
$(MST_DIR):
	mkdir -p $(MST_DIR)

# Rule for generating sequential executable
$(SEQ_EXEC): $(SEQ_OBJS)
	$(CC) $(CFLAGS) $(SEQ_OBJS) -o $(SEQ_EXEC)

# Rule for generating parallel executable
$(PAR_EXEC): $(PAR_OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBDIR) $(PAR_OBJS) -o $(PAR_EXEC) $(LIBS)

# Rule for generating sequential object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule for generating parallel object files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Rule for cleaning files generated
clean:
	rm -f $(SEQ_OBJS) $(PAR_OBJS) $(SEQ_EXEC) $(PAR_EXEC) *.o

# Rule for cleaning all generated mst
clean_mst:
	rm -f ../mst/*.txt

# Rule for cleaning all
clean_all: clean clean_mst