TARGET = halflife
CC = gcc
CFLAGS = -Wall -lm -O2 -pedantic

.PHONY: default all clean

all: $(TARGET)

default: $(TARGET)

$(TARGET): halflife.c
	$(CC) $< $(CFLAGS) -o $@ 
clean:
	\rm -f $(TARGET)

