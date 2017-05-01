CC=gcc
CFLAGS = -O3 -g -Wall -Wextra -pedantic -Wwrite-strings
LDFLAGS = -lz -lm -lpthread -Wall
SOURCES = sonic.c sonic.h sonic_interval.c sonic_interval.h sonic_reference.c sonic_reference.h sonic_structures.h
EXESOURCES = sonic_exe.c
EXEFILE = sonic
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = libsonic.a
INSTALLPATH = /usr/local/lib


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

exe: $(EXESOURCES) $(EXECUTABLE)
	$(CC) $(EXESOURCES) $(EXECUTABLE) $(LDFLAGS) -o $(EXEFILE)

$(EXECUTABLE): $(OBJECTS) 
	ar -rc $(EXECUTABLE)  $(OBJECTS) 

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 

