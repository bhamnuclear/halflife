program1 = halflife
program2 = halflife_double
CC = gcc
CFLAGS = -Wall -lm -O2 -pedantic

.PHONY: default all clean
.DEFAULT_GOAL:=all

all: $(program1) $(program2)

$(program1): %: %.c
	$(CC) $< $(CFLAGS) -o $@ 

$(program2): %: %.c
	$(CC) $< $(CFLAGS) -o $@ 

clean:
	\rm -f $(program1) $(program2)
