CC=mpicc
CFLAGS=-Wall -Wextra
OBJECTS=main.o vp_stack.o vp_master_buffer.o

all: vp

vp: $(OBJECTS)
	$(CC) -o ./vp $(OBJECTS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $^

clean:
	rm -rf *.o vp
