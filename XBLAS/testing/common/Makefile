include ../../make.inc

SRCS = dummy_main.c
OBJS = $(SRCS:.c=.o)

all: $(OBJS)

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o core

