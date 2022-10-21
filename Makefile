CC=gcc
CFLAGS=-Wall -O2

all: synthbar

synthbar: synthbar.c kseq.h
	$(CC) $(CFLAGS) synthbar.c -o $@ -lz

clean:
	rm -rf synthbar *.o
