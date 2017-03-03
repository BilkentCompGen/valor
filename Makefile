CC=gcc
CFLAGS = -O3 -g  
LDFLAGS = -lz -lm -lpthread
SOURCES = sonic.c sonic.h sonic_interval.c sonic_interval.h
TESTSOURCES = testsonic.c
TESTEXE = testsonic
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

