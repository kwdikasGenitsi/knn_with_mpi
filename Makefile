CC=mpicc
CFLAGS=-Wall -Wextra
OBJECTS=main.o

all: vp

vp: $(OBJECTS)
	$(CC) -o ./vp $(OBJECTS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $^

clean:
	rm -rf *.o vp