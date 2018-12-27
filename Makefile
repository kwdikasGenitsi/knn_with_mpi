CC=mpicc
CFLAGS=-Wall -Wextra -g 
OBJECTS=main.o stack.o vp_master_buffer.o array.o vp_tree_local.o point_and_dataset.o

all: vp

vp: $(OBJECTS)
	$(CC) -o ./vp $(OBJECTS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $^

clean:
	rm -rf *.o vp
