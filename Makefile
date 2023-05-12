CC=gcc
CFLAGS=-Wall -O2

all: synthbar

synthbar: synthbar.c kstring.o
	$(CC) $(CFLAGS) $^ -o $@ -lz

kstring.o:
	$(CC) -c $(FLAGS) kstring.c -o $@

clean:
	rm -rf synthbar *.o
