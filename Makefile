VALOR_VERSION := "2.1.5"

VALOR_UPDATE := "05 May 2019"

VALOR_DEBUG := 0
LIVE_PROG := 0
BUILD_DATE := "$(shell date)"
CC=gcc
OPT?= -O2
CFLAGS = -fopenmp  -Wall  $(OPT) -I htslib -I vh -I sonic -DVALOR_VERSION=\"$(VALOR_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DVALOR_UPDATE=\"$(VALOR_UPDATE)\" -DVALOR_DEBUG=$(VALOR_DEBUG) -DLIVE_PROGRESS=$(LIVE_PROG) 
LDFLAGS = -fopenmp htslib/libhts.a -lz -lm -lpthread sonic/libsonic.a
SOURCES =  valor.c cluster.c clique.c clique_inter.c bitset.c hashtable.c statistics.c graph.c common.c vector.c set.c interval10X.c structural_variation.c interc_sv.c cnv.c readbam.c readbed.c recovermolecules.c progress.c cmdline.c config.c valorconfig.c 
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = valor
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

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
