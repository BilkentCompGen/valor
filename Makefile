CC=gcc
CFLAGS = -O0 -g
LDFLAGS = -lz -lm -lpthread
SOURCES = sonic.c sonic.h sonic_interval.c sonic_interval.h sonic_reference.c sonic_reference.h sonic_structures.h
TESTSOURCES = sonic_exe.c
TESTEXE = sonic
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = libsonic.a
INSTALLPATH = /usr/local/lib


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

test: $(TESTSOURCES) $(EXECUTABLE)
	$(CC) $(TESTSOURCES) $(EXECUTABLE) $(LDFLAGS) -o $(TESTEXE)

$(EXECUTABLE): $(OBJECTS) 
	ar -rc $(EXECUTABLE)  $(OBJECTS) 

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 

