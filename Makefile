<<<<<<< HEAD
VALOR_VERSION := "2.1"


VALOR_UPDATE := "24 Sep 2018"


VALOR_DEBUG := 0
LIVE_PROG := 0
BUILD_DATE := "$(shell date)"
CC=gcc
OPT=-O3
CFLAGS = -fopenmp  -Wall  $(OPT) -I htslib -I vh -I sonic -DVALOR_VERSION=\"$(VALOR_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DVALOR_UPDATE=\"$(VALOR_UPDATE)\" -DVALOR_DEBUG=$(VALOR_DEBUG) -DLIVE_PROGRESS=$(LIVE_PROG) 
LDFLAGS = -fopenmp htslib/libhts.a -lz -lm -lpthread sonic/libsonic.a
SOURCES =  valor.c cluster.c clique.c bitset.c hashtable.c statistics.c graph.c common.c vector.c set.c interval10X.c structural_variation.c cnv.c readbam.c readbed.c recovermolecules.c progress.c cmdline.c config.c 
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = valor
INSTALLPATH = /usr/local/bin/
=======
CC=gcc
CFLAGS = -O3 -g -Wall -Wextra -pedantic -Wwrite-strings
LDFLAGS = -lz -lm -lpthread -Wall
SOURCES = sonic.c sonic.h sonic_interval.c sonic_interval.h sonic_reference.c sonic_reference.h sonic_structures.h
EXESOURCES = sonic_exe.c
EXEFILE = sonic
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = libsonic.a
INSTALLPATH = /usr/local/lib
>>>>>>> dab878aca4f5fcae0149d3cd1a62e38eb3187cb3


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

<<<<<<< HEAD
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(OPT) -o $@ $(LDFLAGS)

.c.o:
	@$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~

sonic/sonic.c:
	git clone https://github.com/calkan/sonic.git sonic/
libs: sonic/sonic.c
	make -C htslib
	make -C sonic
install:
	cp $(EXECUTABLE) $(INSTALLPATH)
=======
exe: $(EXESOURCES) $(EXECUTABLE)
	$(CC) $(EXESOURCES) $(EXECUTABLE) $(LDFLAGS) -o $(EXEFILE)

$(EXECUTABLE): $(OBJECTS) 
	ar -rc $(EXECUTABLE)  $(OBJECTS) 

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 
>>>>>>> dab878aca4f5fcae0149d3cd1a62e38eb3187cb3

