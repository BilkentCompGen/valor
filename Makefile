VALOR_VERSION := "0.5-alpha"
VALOR_UPDATE := "16 Mar 2018"
OUT_DIR := "VALOR_OUTPUT"
LOG_FILE := "valor.log"
VALOR_DEBUG := 0
LIVE_PROG := 0
BUILD_DATE := "$(shell date)"
CC=clang
OPT=-O3
CFLAGS =  -Wall  $(OPT) -I htslib -I vh -I sonic -DVALOR_LOG_FILE=\"$(LOG_FILE)\" -DVALOR_VERSION=\"$(VALOR_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DVALOR_UPDATE=\"$(VALOR_UPDATE)\" -DVALOR_DEBUG=$(VALOR_DEBUG) -DLIVE_PROGRESS=$(LIVE_PROG) -DOUT_DIR=\"$(OUT_DIR)\"
LDFLAGS = htslib/libhts.a -lz -lm -lpthread sonic/libsonic.a
SOURCES =  valor.c cluster.c clique.c hashtable.c statistics.c graph.c common.c vector.c set.c interval10X.c structural_variation.c cnv.c readbam.c readbed.c recovermolecules.c progress.c cmdline.c config.c 
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = valor
INSTALLPATH = /usr/local/bin/


all: $(SOURCES) $(EXECUTABLE)
	mkdir -p $(OUT_DIR)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(OPT) -o $@ $(LDFLAGS)

.c.o:
	@$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~

libs:
	make -C htslib
	make -C sonic
install:
	cp $(EXECUTABLE) $(INSTALLPATH)

