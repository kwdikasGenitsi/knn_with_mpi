CC=mpicc
CFLAGS=-Wall -Wextra -Werror -g
OBJECTS=main.o \
	stack.o \
	array.o \
	vp_tree_local.o \
	dataset.o \
	mpi_partition.o \
	vp_tree_distributed.o

all: vp

vp: $(OBJECTS)
	$(CC) -o ./vp $(OBJECTS) -lm

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $^

clean:
	rm -rf *.o vp
