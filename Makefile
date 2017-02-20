CC=gcc
CFLAGS = -O3 -g  
LDFLAGS = -lz -lm -lpthread
SOURCES = sonic.c sonic.h 
	
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = libsonic.a
INSTALLPATH = /usr/local/lib

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	ar -rc $(EXECUTABLE)  $(OBJECTS) 

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 

